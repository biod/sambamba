/*
    This file is part of Sambamba.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
/** module for executing samtools mpileup in parallel using named pipes,
 *  after chunking a file 
 */

module sambamba.pileup;

import sambamba.utils.common.bed;
import sambamba.utils.common.tmpdir;

import bio.bam.multireader;
import bio.bam.reader;
import bio.bam.writer;
import bio.bam.pileup;

import bio.core.utils.format : write;
import bio.core.utils.roundbuf;

import std.process;
import std.stdio;
import std.parallelism;
import std.file : rmdirRecurse;
import std.algorithm;
import std.array;
import std.getopt;
import std.string : strip, indexOf, toStringz;
import std.c.stdlib;
import std.typecons;
import std.stream;
import std.range;
import std.algorithm;
import std.path;
import std.traits;
import std.typecons;
import std.conv;

import core.thread;
import core.sync.mutex;
import core.sys.posix.sys.stat;
import core.sys.posix.stdio : fopen;
import core.stdc.errno;

extern(C) char* mkdtemp(char* template_);
extern(C) int mkfifo(immutable(char)* fn, int mode);

void makeFifo(string filename) {
    auto s = toStringz(filename);
    int ret = mkfifo(s, octal!"666");
    if (ret == -1) {
        stderr.writeln(errno);
        throw new Exception("failed to create named pipe " ~ filename);
    }
}

struct BedRecord {
    string reference;
    ulong start;
    ulong end;

    void toString(scope void delegate(const(char)[]) dg) const {
        dg.write(reference);
        dg.write('\t');
        dg.write(start);
        dg.write('\t');
        dg.write(end);
    }
}

struct ForkData(C) {
    MultiBamReader bam;
    TaskPool task_pool;
    C chunk;
    string filename;

    void toString(scope void delegate(const(char)[]) dg) const {
        dg.write(filename);
        dg.write('\n');
        dg.write(bam.reference_sequences[chunk.ref_id].name);
        dg.write('\n');
        dg.write(chunk.start_position);
        dg.write('\n');
        dg.write(chunk.end_position);
    }

    const(BedRecord) bed() @property const {
        return BedRecord(bam.reference_sequences[chunk.ref_id].name,
                         chunk.start_position + 1,
                         chunk.end_position);
    }
}

struct MArray(T) { T[] data; T* ptr; }

// TODO: detect if samtools writes in pileup format
MArray!char runSamtools(string filename,
                        string[] samtools_args, string[] bcftools_args)
{
    auto samtools_cmd = (["samtools mpileup", filename, "-gu",
                          "-l", filename ~ ".bed"] ~ samtools_args).join(" ");
    auto bcftools_cmd = (["bcftools view -"] ~ bcftools_args).join(" ");
    auto cmd = samtools_cmd ~ " | " ~ bcftools_cmd;

    stderr.writeln("[executing] ", cmd);
    auto pp = pipeShell(cmd, Redirect.stdout | Redirect.stderr);
    scope(exit) pp.pid.wait();

    size_t capa = 1_024_576;
    size_t used = 0;
    char* result = cast(char*)std.c.stdlib.malloc(capa);

    char[4096] buffer = void;
    while (true) {
        auto buf = pp.stdout.rawRead(buffer[]);
        if (buf.length == 0)
            break;
        if (used + buf.length > capa) {
            capa = max(capa * 2, used + buf.length);
            result = cast(char*)std.c.stdlib.realloc(cast(void*)result, capa);
        }
        result[used .. used + buf.length] = buf[];
        used += buf.length;
    }

    auto vcf = result[0 .. used];
    if (!filename.endsWith("/1")) { 
        // if it's not the first file
        // strip lines starting with '#'
        while (true) {
            auto pos = vcf.indexOf('\n');
            if (pos != -1 && pos + 1 < vcf.length && vcf[pos + 1] == '#')
                vcf = vcf[pos + 1 .. $];
            else
                break;
        }

        if (vcf.startsWith("#")) {
            auto pos = vcf.indexOf('\n');
            if (pos != -1)
                vcf = vcf[pos + 1.. $];
            else
                vcf = vcf[$ .. $];
        }
    }

    return typeof(return)(vcf, result);
}

class PileupRunner {
    alias Task!(runSamtools, string, string[], string[])* PileupTask;
    TaskPool task_pool;
    PileupTask[] tasks;
    string[] samtools_args;
    string[] bcftools_args;
    std.stdio.File output_file;

    size_t finished, submitted;

    this(TaskPool pool, std.stdio.File output_file) {
        task_pool = pool;
        this.output_file = output_file;
        tasks = new PileupTask[task_pool.size * 2];
    }

    void dumpNextVcf() {
        auto vcf = tasks[finished % $].yieldForce();
        tasks[finished++ % $] = null;
        output_file.rawWrite(vcf.data);
        std.c.stdlib.free(vcf.ptr);
    }
    
    void pushTask(string filename) {
        if (submitted - finished == tasks.length)
            dumpNextVcf();

        tasks[submitted % $] = task!runSamtools(filename,
                                                samtools_args, bcftools_args);
        task_pool.put(tasks[submitted++ % $]);

        while (finished < submitted && tasks[finished % $].done()) {
            dumpNextVcf();
        }
    }

    void wait() {
        while (finished != submitted)
            dumpNextVcf();
    }
}

class ChunkDispatcher(ChunkRange) {
    private string tmp_dir_;
    private ChunkRange chunks_;
    private MultiBamReader bam_;
    private size_t num_;
    private Mutex mutex_;
    private PileupRunner runner_;

    alias ElementType!(Unqual!(ChunkRange)) Chunk;

    this(string tmp_dir, ChunkRange chunks, MultiBamReader bam,
         PileupRunner runner) 
    {
        tmp_dir_ = tmp_dir;
        chunks_ = chunks;
        bam_ = bam;
        num_ = 0;
        mutex_ = new Mutex();
        runner_ = runner;
    }

    Nullable!(Tuple!(Chunk, string)) nextChunk() {
        mutex_.lock();
        scope(exit) mutex_.unlock();

        typeof(return) chunk;
        if (chunks_.empty)
            return chunk;
        ++num_;

        auto filename = buildPath(tmp_dir_, num_.to!string()); 
        makeFifo(filename);
        chunk = tuple(chunks_.front, filename);
        chunks_.popFront();

        auto ref_name = bam_.reference_sequences[chunk[0].ref_id].name;
        auto start = chunk[0].start_position + 1;
        auto end = chunk[0].end_position;
            
        auto bed = BedRecord(ref_name, start - 1, end);
        std.stdio.File(filename ~ ".bed", "w").writeln(bed);

        runner_.pushTask(filename);

        return chunk;
    }
}

void worker(Dispatcher)(Dispatcher d,
                        MultiBamReader bam,
                        TaskPool task_pool) {
    while (true) {
        auto result = d.nextChunk();
        if (result.isNull)
            return;

        auto chunk = result[0];
        auto filename = result[1];

        auto writer = new BamWriter(filename, 0, task_pool);
        scope(exit) writer.finish();
        with (writer) {
            writer.writeSamHeader(bam.header);
            writer.writeReferenceSequenceInfo(bam.reference_sequences);
            foreach (read; chunk.reads)
                writer.writeRecord(read);
        }
    }
}

auto chunkDispatcher(ChunkRange)(string tmp_dir, ChunkRange chunks, 
                                 MultiBamReader bam, PileupRunner runner) {
    return new ChunkDispatcher!ChunkRange(tmp_dir, chunks, bam, runner);
}

void printUsage() {
    stderr.writeln("usage: sambamba-pileup [options] input.bam [input2.bam [...]]");
    stderr.writeln("                       [--samtools <samtools mpileup args>]");
    stderr.writeln("                       [--bcftools <bcftools args>]");
    stderr.writeln();
    stderr.writeln("This subcommand relies on external tools and acts ");
    stderr.writeln("as a wrapper boosting the performance by parallelization,");
    stderr.writeln("rather than trying to duplicate or introduce new algorithms");
    stderr.writeln("for variant calling.");
    stderr.writeln("The following tools should be present in $PATH:");
    stderr.writeln("    * samtools");
    stderr.writeln("    * bcftools");
    stderr.writeln();
    stderr.writeln("If --bcftools is skipped, bcftools is not called");
    stderr.writeln("If --samtools is skipped, samtools mpileup is called with default arguments -gu -S -D -d 1000 -L 1000 -m 3 -F 0.0002");
    stderr.writeln();
    stderr.writeln("Splits input into multiple chunks and feeds them");
    stderr.writeln("to samtools mpileup and, optionally, bcftools in parallel.");
    stderr.writeln("The chunks are slightly overlapping so that variant calling");
    stderr.writeln("should not be impacted by these manipulations. The obtained results");
    stderr.writeln("from the multiple processes are then combined into a single file.");
    stderr.writeln("Options: -F, --filter=FILTER");
    stderr.writeln("                    set custom filter for alignments");
    stderr.writeln("         -L, --regions=FILENAME");
    stderr.writeln("                    provide BED file with regions");
    stderr.writeln("                    (no need to duplicate it in samtools args);");
    stderr.writeln("                    all input files must be indexed");
    stderr.writeln("         -o, --output-filename=<STDOUT>");
    stderr.writeln("                    specify output filename");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -b, --buffer-size=4_000_000");
    stderr.writeln("                    chunk size (in bytes)");
}

version(standalone) {
    int main(string[] args) {
       return pileup_main(args);
    }
}

int pileup_main(string[] args) {
    auto bcftools_args = find(args, "--bcftools");
    auto args1 = (bcftools_args.length>0 ? args[0 .. $-bcftools_args.length] : args );
    auto samtools_args = find(args1, "--samtools");
    auto args2 = (samtools_args.length>0 ? args1[0 .. $-samtools_args.length] : args1 );

    if (!samtools_args.empty) {
        samtools_args.popFront();
    } else {
        // Default values for samtools if not passed 
        samtools_args = ["-gu",          // uncompressed BCF output
                         "-S",           // per-sample strand bias
                         "-D",           // per-sample DP
                         "-d", "1000",   // max per-BAM depth
                         "-L", "1000",   // max depth for indel calling
                         "-m", "3",      // min. gapped reads for indel candidates
                         "-F", "0.0002", // min. fraction of gapped reads
                         ];
    }

    if (!bcftools_args.empty) {
        bcftools_args.popFront();
    }

    auto own_args = args2;

    string bed_filename;
    string query;
    string output_filename;
    uint n_threads = defaultPoolThreads;
    std.stdio.File output_file = stdout;
    size_t buffer_size = 4_000_000;

    try {
        getopt(own_args,
               std.getopt.config.caseSensitive,
               "regions|L",         &bed_filename,
               "filter|F",          &query,
               "output-filename|o", &output_filename,
               "nthreads|t",        &n_threads,
               "buffer-size|b",     &buffer_size);

        if (own_args.length < 2) {
            printUsage();
            return 0;
        }

        defaultPoolThreads = n_threads;
        auto bam = new MultiBamReader(args[1 .. $]);

        char[] buf = defaultTmpDir() ~ "/sambamba-fork-XXXXXX\0".dup;
        mkdtemp(buf.ptr);

        string tmp_dir = to!string(buf.ptr);
        scope(exit) rmdirRecurse(tmp_dir);

        auto pileupRunner = new PileupRunner(taskPool, output_file);
        pileupRunner.samtools_args = samtools_args;
        pileupRunner.bcftools_args = bcftools_args;

        InputRange!BamRead reads;
        if (bed_filename is null) {
            reads = inputRangeObject(bam.reads().map!`a.read`);
        } else {
            auto regions = parseBed(bed_filename, bam);
            reads = inputRangeObject(bam.getReadsOverlapping(regions).map!`a.read`);
        }

        auto chunks = bam.reads().pileupChunks(false, buffer_size);
        auto dispatcher = chunkDispatcher(tmp_dir, chunks, bam, pileupRunner);
    
        auto threads = new ThreadGroup();
        foreach (i; 0 .. n_threads)
            threads.create(() { worker(dispatcher, bam, taskPool); });

        threads.joinAll();
        pileupRunner.wait();
        return 0;

    } catch (Exception e) {
        stderr.writeln("sambamba-pileup: ", e.msg);

        version(development) {
            throw e;
        }

        return 1;
    }

    return 0;
}
