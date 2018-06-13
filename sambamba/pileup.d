/*
    This file is part of Sambamba.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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
import sambamba.utils.common.overwrite;
import utils.lz4;

import bio.bam.multireader;
import bio.bam.reader;
import bio.bam.writer;
import bio.bam.pileup;

import bio.core.utils.format : write;
import bio.core.utils.roundbuf;
import bio.core.utils.stream;

import std.process;
import std.stdio;
import std.parallelism;
import std.file : rmdirRecurse, exists;
import std.algorithm;
import std.array;
import std.exception;
import std.getopt;
import std.string : strip, indexOf, toStringz;
import core.stdc.stdlib;
import std.typecons;
import undead.stream;
import std.range;
import std.algorithm;
import std.path;
import std.regex;
import std.traits;
import std.typecons;
import std.conv;
import std.container;

import core.thread;
import core.sync.mutex, core.sync.condition;
import core.sys.posix.sys.stat;
import core.sys.posix.stdio : fopen;
import core.stdc.errno;

extern(C) char* mkdtemp(char* template_);
extern(C) int mkfifo(immutable(char)* fn, int mode);

// Cached values
string samtoolsBin     = null;
string samtoolsVersion = null;
string bcftoolsBin     = null;
string bcftoolsVersion = null;

// Return path to samtools after testing whether it exists and supports mpileup
auto samtoolsInfo()
{
  if (samtoolsBin is null) {
    auto paths = environment["PATH"].split(":");
    auto a = array(filter!(path => exists(path ~ "/samtools"))(paths));
    enforce(!a.empty, "failed to locate samtools executable in PATH");
    samtoolsBin = a[0] ~ "/samtools";
  }
  enforce(exists(samtoolsBin), samtoolsBin ~ " is invalid");
  if (samtoolsVersion is null) {
    auto samtools = execute([samtoolsBin]);
    enforce(samtools.status==1, "samtools failed: " ~ samtools.output);
    auto r = regex(r"Version: 1\.\d(\.\d)?[^\n]+");
    enforce(matchFirst(samtools.output, r), "Can not find version in " ~ samtools.output);
    samtoolsVersion = matchFirst(samtools.output, r).hit;
    enforce(samtoolsVersion.startsWith("Version: 1."), "version " ~ samtoolsVersion ~ " of samtools is unsupported");
  }
  return [samtoolsBin, samtoolsVersion];
}

auto samtoolsPath() { return samtoolsInfo()[0]; }

auto bcftoolsInfo()
{
  if (bcftoolsBin is null) {
    auto paths = environment["PATH"].split(":");
    auto a = array(filter!(path => exists(path ~ "/bcftools"))(paths));
    enforce(!a.empty, "failed to locate bcftools executable in PATH");
    bcftoolsBin = a[0] ~ "/bcftools";
  }
  enforce(exists(bcftoolsBin), bcftoolsBin ~ " is invalid");
  if (bcftoolsVersion is null) {
    auto bcftools = execute([bcftoolsBin]);
    enforce(bcftools.status == 1, "bcftools failed: " ~ bcftools.output);
    auto r = regex(r"Version: 1\.\d(\.\d)?[^\n]+");
    enforce(matchFirst(bcftools.output, r), "Can not find version in " ~ bcftools.output);
    bcftoolsVersion = matchFirst(bcftools.output, r).hit;
    enforce(bcftoolsVersion.startsWith("Version: 1."), "version " ~ bcftoolsVersion ~ " of bcftools is unsupported");
  }
  return [bcftoolsBin, bcftoolsVersion];
}

auto bcftoolsPath() { return bcftoolsInfo()[0]; }

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

struct MArray(T) { T[] data; T* ptr; }

MArray!char data;

__gshared string this_app;

private {
    struct Recipe {
        string strip_header_cmd;
        string compression_cmd;
        void function (ubyte[] data, std.stdio.File output_file) decompressionFunc;
    }

    void dump(ubyte[] data, std.stdio.File output_file) {
        output_file.rawWrite(data);
    }

    void lz4decompress(ubyte[] data, std.stdio.File output_file) {
        lz4decompressor.decompress(new MemoryStream(data), output_file);
    }

    __gshared Recipe[FileFormat] recipes;
    __gshared LZ4Decompressor lz4decompressor;
}

// TODO: fix bcftoolsCommand and samtoolsCommand in the header
void init() {
    lz4decompressor = new LZ4Decompressor();

    recipes[FileFormat.pileup] =          Recipe(this_app ~ " strip_bcf_header --vcf",
                                                 this_app ~ " lz4compress",
                                                 &lz4decompress);
    recipes[FileFormat.BCF] =             Recipe(this_app ~ " strip_bcf_header --bcf",
                                                 null,
                                                 &dump);
    recipes[FileFormat.uncompressedBCF] = Recipe(this_app ~ " strip_bcf_header --ubcf",
                                                 this_app ~ " lz4compress",
                                                 &lz4decompress);
    recipes[FileFormat.VCF] =             Recipe(this_app ~ " strip_bcf_header --vcf",
                                                 this_app ~ " lz4compress",
                                                 &lz4decompress);
}

string makeInputCmdLine(string input_cmd, FileFormat input_format, bool strip_header) {
    auto recipe = recipes[input_format];
    string cmd = input_cmd;
    if (strip_header && recipe.strip_header_cmd !is null)
        cmd ~= "| " ~ recipe.strip_header_cmd;
    if (recipe.compression_cmd !is null)
        cmd ~= "| " ~ recipe.compression_cmd;
    return cmd;
}

void decompressIntoFile(char[] data, FileFormat input_format,
                        std.stdio.File output_file) {
    recipes[input_format].decompressionFunc(cast(ubyte[])data, output_file);
}

struct Args {
    string[] samtools_args;
    string[] bcftools_args;
    FileFormat input_format;

    this(string[] samtools_args_, bool use_bcftools, string[] bcftools_args_) {
        samtools_args = unbundle(samtools_args_);
        bcftools_args = unbundle(bcftools_args_, "O"); // keep -Ov|-Ob|...
        auto samtools_output_fmt = fixSamtoolsArgs(use_bcftools, samtools_args);
        auto bcftools_output_fmt = fixBcftoolsArgs(use_bcftools, bcftools_args, samtools_args);

        input_format = samtools_output_fmt;
        if (use_bcftools)
            input_format = bcftools_output_fmt;
    }

    string makeCommandLine(string filename) {
        auto basic_args = [samtoolsPath(), "mpileup", filename];
        basic_args ~= ["-l", filename ~ ".bed"];
        auto samtools_cmd = (basic_args ~ samtools_args).join(" ");
        string cmd = samtools_cmd;
        if (bcftools_args.length > 0) {
            auto bcftools_cmd = bcftoolsPath() ~ " " ~ bcftools_args.join(" ");
            cmd = samtools_cmd ~ " | " ~ bcftools_cmd;
        }

        bool strip_header = !filename.endsWith("/1");
        return makeInputCmdLine(cmd, input_format, strip_header);
    }
}

enum FileFormat {
    pileup,
    BCF,
    uncompressedBCF,
    VCF,
    gzippedVCF
}

string[] unbundle(string[] args, string exclude="") {
    import std.ascii : isAlpha;
    import std.conv : text;
    import std.algorithm : count;
    string[] unbundled;
    foreach (a; args) {
        if (a.length >= 2 && a[0] == '-' && exclude.count(a[1]) == 0) {
            string[] expanded;
            foreach (j, dchar c; a[1 .. $])
            {
                if (!isAlpha(c)) {
                    expanded ~= a[j + 1 .. $];
                    break;
                }
                expanded ~= text('-', c);
            }
            unbundled ~= expanded;
        } else {
            unbundled ~= a;
        }
    }
    return unbundled;
}

// input: unbundled samtools arguments
// output: detected output format
FileFormat fixSamtoolsArgs(bool use_bcftools, ref string[] args) {
    bool vcf = false;
    bool bcf = false;
    bool uncompressed = false;
    bool[] keep;
    foreach (i; 0 .. args.length) {
        if (args[i] == "-o") {
            throw new Exception("-o argument of samtools is disallowed, use --output-filename argument of sambamba mpileup");
        }
        if (args[i] == "-g") {
            bcf = true; keep ~= true;
        } else if (args[i] == "-v") {
            vcf = true; keep ~= !use_bcftools;
        } else if (args[i] == "-u") {
            bcf = true; uncompressed = true; keep ~= true;
        } else {
            keep ~= true;
        }
    }

    string[] fixed_args;
    foreach (i; 0 .. args.length) {
        if (keep[i])
            fixed_args ~= args[i];
    }

    // When using bcftools add these switches
    bool fixes_applied;
    if (vcf && use_bcftools) {
        fixed_args ~= ["-g", "-u"];
        fixes_applied = true;
    } else if (bcf && use_bcftools && !uncompressed) {
        fixed_args ~= "-u";
        fixes_applied = true;
    }

    args = fixed_args;

    if (fixes_applied && use_bcftools) {
        stderr.writeln("NOTE: changed samtools output format to uncompressed BCF for better performance (-gu)");
    }

    if (bcf && vcf) {
        throw new Exception("samtools can't be asked for both -g and -v");
    } else if (bcf && uncompressed) {
        return FileFormat.uncompressedBCF;
    } else if (bcf && !uncompressed) {
        return FileFormat.BCF;
    } else if (vcf && uncompressed) {
        return FileFormat.VCF;
    } else if (vcf && !uncompressed) {
        // TODO
        throw new Exception("compressed VCF is not supported, please use bgzip and uncompressed VCF");
    } else {
        return FileFormat.pileup;
    }
}

// input: unbundled bcftools arguments
// output: detected output format
FileFormat fixBcftoolsArgs(bool use_bcftools, ref string[] args, ref string[] samtools_args) {
    FileFormat fmt = FileFormat.VCF;
    bool[] keep;
    foreach (i; 0 .. args.length) {
        if (args[i] == "-o") {
            throw new Exception("-o argument of bcftools is disallowed, use --output-filename argument of sambamba mpileup");
        }
        if (args[i] == "-Ov") {
            fmt = FileFormat.VCF; keep ~= true;
        } else if (args[i] == "-Oz") {
            // TODO
            throw new Exception("compressed VCF is not supported, please use bgzip and uncompressed VCF");
            // fmt = FileFormat.gzippedVCF; keep ~= false;
        } else if (args[i] == "-Ob") {
            fmt = FileFormat.BCF; keep ~= true;
        } else if (args[i] == "-Ou") {
            fmt = FileFormat.uncompressedBCF; keep ~= true;
        } else {
            keep ~= true;
        }
    }

    string[] fixed_args;
    foreach (i; 0 .. args.length) {
        if (keep[i])
            fixed_args ~= args[i];
    }
    // When using bcftools and args is empty add these switches
    if (use_bcftools && fixed_args.empty) {
      fixed_args = ["view", "-"];
      samtools_args ~= [ "-g", "-u" ];
    }
    args = fixed_args;
    return fmt;
}

class ChunkDispatcher(ChunkRange) {
    private string tmp_dir_;
    private ChunkRange chunks_;
    private MultiBamReader bam_;
    private size_t num_, curr_num_ = 1;
    private int total_num_;
    private FileFormat format_;
    private std.stdio.File output_file_;
    private size_t max_queue_length_;
    private int prev_ref_id;
    private ulong prev_pos_diff;
    private size_t max_queue_data_size_;
    private size_t curr_queue_data_size_;

    private Mutex mutex_, queue_mutex_;
    private Condition queue_not_empty_condition_, queue_not_full_condition_;

    alias Tuple!(size_t, "num", char[], "data") Result;
    alias Array!(Result) ResultQueue;
    private BinaryHeap!(ResultQueue, "a.num > b.num") result_queue_;

    alias ElementType!(Unqual!(ChunkRange)) Chunk;

    this(string tmp_dir, ChunkRange chunks, MultiBamReader bam,
         FileFormat format, std.stdio.File output_file, size_t max_queue_length,
         size_t max_queue_data_size) {
        tmp_dir_ = tmp_dir;
        chunks_ = chunks;
        bam_ = bam;
        num_ = 0;
        total_num_ = -1;
        mutex_ = new Mutex();
        queue_mutex_ = new Mutex();
        queue_not_empty_condition_ = new Condition(queue_mutex_);
        queue_not_full_condition_ = new Condition(queue_mutex_);
        result_queue_ = heapify!("a.num > b.num", ResultQueue)(ResultQueue());
        format_ = format;
        output_file_ = output_file;
        max_queue_length_ = max_queue_length;
        max_queue_data_size_ = max_queue_data_size;
    }

    Nullable!(Tuple!(Chunk, string, size_t)) nextChunk() {
        mutex_.lock();
        scope(exit) mutex_.unlock();

        typeof(return) chunk;
        if (chunks_.empty) {
            if (num_ == 0)
                total_num_ = num_.to!int();
            return chunk;
        }

        if (chunks_.front.ref_id >= bam_.reference_sequences.length) {
            if (chunks_.front.ref_id != -1)
              stderr.writeln("Invalid ref. id ", chunks_.front.ref_id, " in chunk #", num_ + 1);
            return chunk;
        }
        ++num_;

        auto filename = buildPath(tmp_dir_, num_.to!string());
        chunk = tuple(chunks_.front, filename, num_);
        chunks_.popFront();
        if (chunks_.empty) {
            stderr.writeln("[Last chunk fetched] ", num_);
            total_num_ = num_.to!int();
        }

        if (num_ > 1) {
            ulong diff = chunk[0].end_position - chunk[0].start_position;
            int ref_id = chunk[0].ref_id;
            if (ref_id == prev_ref_id && !chunk[0].front.reads.empty &&
                prev_pos_diff + diff < chunk[0].front.reads[0].sequence_length) // assume reads are ~same length
                stderr.writeln("[WARNING] COVERAGE IS TOO HIGH, INCREASE --buffer-size TO AVOID WRONG RESULTS");
            prev_pos_diff = diff;
            prev_ref_id = ref_id;
        }

        auto ref_name = bam_.reference_sequences[chunk[0].ref_id].name;
        auto f = std.stdio.File(filename ~ ".bed", "w");
        if (bed_filename is null) {
            auto start = chunk[0].start_position;
            auto end = chunk[0].end_position;

            auto bed = BedRecord(ref_name, start, end);
            f.writeln(bed);
        } else {
            foreach (reg; regions) {
                if (chunk[0].ref_id != reg.ref_id) continue;
                auto start = max(reg.start, chunk[0].start_position);
                auto end = min(reg.end, chunk[0].end_position);
                if (start > end) continue;
                auto bed = BedRecord(ref_name, start, end);
                f.writeln(bed);
            }
        }
        f.close();

        return chunk;
    }

    void queueResult(size_t num, char[] data) {
        synchronized(queue_mutex_) {
            while (true) { // wait until it's good time to put the result into the output queue
                if (num == curr_num_) // always put if it should be written out now
                    break;
                if (result_queue_.empty) // always put if the queue is empty
                    break;

                if (result_queue_.length >= max_queue_length_) {
                    stderr.writeln("[chunk waiting for dump queue] ", num, " (output is too slow: reduce threads or improve output speed)");
                    queue_not_full_condition_.wait();
                } else if (curr_queue_data_size_ + data.length > max_queue_data_size_) {
                    stderr.writeln("[chunk waiting for dump queue] ", num, " (output is too large: increase buffer size)");
                    queue_not_full_condition_.wait();
                } else {
                    break;
                }
            }
            result_queue_.insert(Result(num, data));
            curr_queue_data_size_ += data.length;
            queue_not_empty_condition_.notify();
        }
        stderr.writeln("[chunk queued for dumping] ", num);
    }

    void dumpResults() {
        Result result;
        while (!dumpFinished()) {
            synchronized(queue_mutex_) {
                while (!hasDumpableResult()) {
                    queue_not_empty_condition_.wait();
                }
                result = result_queue_.front;
                result_queue_.popFront();
                curr_queue_data_size_ -= result.data.length;
                ++curr_num_;
                queue_not_full_condition_.notifyAll();
            }
            decompressIntoFile(result.data, format_, output_file_);
            stderr.writeln("[chunk dumped] ", result.num);
            core.stdc.stdlib.free(result.data.ptr);
        }
    }

private:
    bool hasDumpableResult() {
        synchronized(queue_mutex_) {
            return !result_queue_.empty() && result_queue_.front.num == curr_num_;
        }
    }

    bool dumpFinished() {
        return total_num_ >= 0 && curr_num_.to!int() > total_num_;
    }
}

void worker(Dispatcher)(Dispatcher d,
                        MultiBamReader bam,
                        Args args) {
    while (true) {
        auto result = d.nextChunk();
        if (result.isNull)
            return;

        auto chunk = result[0];
        auto filename = result[1];
        auto num = result[2];
        makeFifo(filename);

        import core.sys.posix.signal;
        signal(SIGPIPE, SIG_IGN);

        auto writing_thread = new Thread(() {

            auto output_stream = new bio.core.utils.stream.File(filename, "w");
            stderr.writeln("[opened FIFO for writing] ", filename);
            auto writer = new BamWriter(output_stream, 0);
            writer.writeSamHeader(bam.header);
            writer.writeReferenceSequenceInfo(bam.reference_sequences);
            foreach (read; chunk.reads)
                writer.writeRecord(read);
            writer.finish();
            stderr.writeln("[closed FIFO] ", filename);
            });

        auto cmd = args.makeCommandLine(filename);
        stderr.writeln("[executing] ", cmd);
        auto pp = pipeShell(cmd, Redirect.stdout);

        writing_thread.start();

        size_t capa = 1_024_576;
        size_t used = 0;
        char* output = cast(char*)core.stdc.stdlib.malloc(capa);

        char[4096] buffer = void;
        while (true) {
            auto buf = pp.stdout.rawRead(buffer[]);
            if (buf.length == 0)
                break;
            if (used + buf.length > capa) {
                capa = max(capa * 2, used + buf.length);
                output = cast(char*)core.stdc.stdlib.realloc(cast(void*)output, capa);
                if (output is null)
                    throw new Exception("failed to allocate " ~ capa.to!string ~ " bytes");
            }
            output[used .. used + buf.length] = buf[];
            used += buf.length;
        }

        writing_thread.join();
        pp.pid.wait();

        d.queueResult(num, output[0 .. used]);
    }
}

auto chunkDispatcher(ChunkRange)(string tmp_dir, ChunkRange chunks,
                                 MultiBamReader bam, FileFormat format,
                                 std.stdio.File output_file, size_t max_queue_length, size_t max_queue_data_size) {
  return new ChunkDispatcher!ChunkRange(tmp_dir, chunks, bam, format, output_file, max_queue_length, max_queue_data_size);
}

void printUsage() {
    stderr.writeln("usage: sambamba-pileup [options] input.bam [input2.bam [...]]");
    stderr.writeln("                       [--samtools <samtools mpileup args>]");
    stderr.writeln("                       [--bcftools <bcftools call args>]");
    stderr.writeln();
    stderr.writeln("This subcommand relies on external tools and acts as a multi-core");
    stderr.writeln("implementation of samtools and bcftools.");
    stderr.writeln();
    stderr.writeln("Therefore, the following tools should be present in $PATH:");
    stderr.writeln("    * samtools");
    stderr.writeln("    * bcftools (when used)");
    stderr.writeln();
    stderr.writeln("If --samtools is skipped, samtools mpileup is called with default arguments");
    stderr.writeln("If --bcftools is used without parameters, samtools is called with");
    stderr.writeln("     switch '-gu' and bcftools is called as 'bcftools view -'");
    stderr.writeln("If --bcftools is skipped, bcftools is not called");
    stderr.writeln();
    stderr.writeln("Sambamba splits input BAM files into chunks and feeds them");
    stderr.writeln("to samtools mpileup and, optionally, bcftools in parallel.");
    stderr.writeln("The chunks are slightly overlapping so that variant calling");
    stderr.writeln("should not be impacted by these manipulations. The obtained results");
    stderr.writeln("from the multiple processes are combined as ordered output.");
    stderr.writeln();
    stderr.writeln("Sambamba options:");
//    stderr.writeln("         -F, --filter=FILTER");
//    stderr.writeln("                    set custom filter for alignments");
    stderr.writeln("         -L, --regions=FILENAME");
    stderr.writeln("                    provide BED file with regions");
    stderr.writeln("                    (no need to duplicate it in samtools args);");
    stderr.writeln("                    all input files must be indexed");
    stderr.writeln("         -o, --output-filename=<STDOUT>");
    stderr.writeln("                    specify output filename");
    stderr.writeln("         --tmpdir=TMPDIR");
    stderr.writeln("                    directory for temporary files");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -b, --buffer-size=64_000_000");
    stderr.writeln("                    chunk size (in bytes)");
    stderr.writeln("         -B, --output-buffer-size=512_000_000");
    stderr.writeln("                    output buffer size (in bytes)");
    stderr.writeln();
    stderr.writeln("Sambamba paths:\n");
    samtoolsInfo();
    stderr.writeln("         samtools: ", samtoolsBin, " ", samtoolsVersion);
    bcftoolsInfo();
    stderr.writeln("         bcftools: ", bcftoolsBin, " ", bcftoolsVersion);
}

version(standalone) {
    int main(string[] args) {
       return pileup_main(args);
    }
}

string output_filename = null;
__gshared string bed_filename = null;
__gshared BamRegion[] regions;

int pileup_main(string[] args) {
    this_app = args[0];
    init();

    auto bcftools_args = find(args, "--bcftools");
    immutable use_bcftools = !bcftools_args.empty;
    auto args1 = (use_bcftools ? args[0 .. $-bcftools_args.length] : args );
    auto samtools_args = find(args1, "--samtools");
    auto own_args = (samtools_args.length>0 ? args1[0 .. $-samtools_args.length] : args1 );

    if (!samtools_args.empty) {
        samtools_args.popFront(); // remove --samtools switch
    } else {
        samtools_args = [];
    }

    if (use_bcftools) {
        bcftools_args.popFront(); // remove the switch --bcftools
    }

    uint n_threads = defaultPoolThreads;
    std.stdio.File output_file = stdout;
    size_t buffer_size = 64_000_000;
    size_t output_buffer_size = 512_000_000;

    string tmp_dir_prefix = defaultTmpDir();

    try {
        getopt(own_args,
               std.getopt.config.caseSensitive,
               "regions|L",         &bed_filename,
               //"filter|F",          &query,
               "output-filename|o", &output_filename,
               "tmpdir",            &tmp_dir_prefix,
               "nthreads|t",        &n_threads,
               "buffer-size|b",     &buffer_size,
               "output-buffer-size|B", &output_buffer_size);

        if (own_args.length < 2) {
            printUsage();
            return 0;
        }

        samtoolsInfo();      // initialize samtools path before threading
        if (use_bcftools)
            bcftoolsInfo();  // initialize bcftools path before threading

        stderr.writeln("samtools mpileup options: ", samtools_args.join(" "));
        if (use_bcftools)
            stderr.writeln("bcftools options: ", bcftools_args.join(" "));

        if (output_filename != null) {
            foreach (filename; own_args[1 .. $])
                protectFromOverwrite(filename, output_filename);
            output_file = std.stdio.File(output_filename, "w+");
        }

        defaultPoolThreads = n_threads;
        auto bam = new MultiBamReader(own_args[1 .. $]);

        string tmp_dir = randomSubdir(tmp_dir_prefix);
        scope(exit) rmdirRecurse(tmp_dir);

        auto bundled_args = Args(samtools_args, use_bcftools, bcftools_args);

        InputRange!BamRead reads;
        if (bed_filename is null) {
            reads = inputRangeObject(bam.reads().map!`a.read`);
        } else {
            regions = parseBed(bed_filename, bam);
            reads = inputRangeObject(bam.getReadsOverlapping(regions).map!`a.read`);
        }

        auto chunks = reads.pileupChunks(false, buffer_size);
        auto dispatcher = chunkDispatcher(tmp_dir, chunks, bam, bundled_args.input_format, output_file,
                                          2 * n_threads, output_buffer_size);

        auto writer = new Thread(&dispatcher.dumpResults);
        writer.start();
        auto threads = new ThreadGroup();

        scope (exit) {
            output_file.close();
        }

        foreach (i; 0 .. max(1, n_threads))
            threads.create(() { worker(dispatcher, bam, bundled_args); });

        threads.joinAll();
        writer.join();
        stderr.writeln("[Successful exit]");

        return 0;

    } catch (Exception e) {
        stderr.writeln("sambamba-pileup: ", e.msg);

        debug {
          throw e;
        }
        else {
          version(development) {
            throw e;
          }
          return 1;
        }
    }
}
