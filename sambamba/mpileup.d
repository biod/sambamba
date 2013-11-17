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
    BamReader bam;
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

MArray!char runSamtools(string filename, string ref_fn) {
    auto cmd = "samtools mpileup " ~ filename ~ " -gu"
               ~ " -l " ~ filename ~ ".bed "
               ~ " -SD -d 1000 -L 1000 -m 3 -F 0.0002 -f " ~ ref_fn
               ~ "| bcftools view -"
//               ~ " -vcg"
               ;
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
    alias Task!(runSamtools, string, string)* PileupTask;
    TaskPool task_pool;
    PileupTask[] tasks;
    string ref_fn;
    std.stdio.File output_file;

    size_t finished, submitted;

    this(TaskPool pool, std.stdio.File output_file, string ref_fn) {
        task_pool = pool;
        this.ref_fn = ref_fn;
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

        tasks[submitted % $] = task!runSamtools(filename, ref_fn);
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
    private BamReader bam_;
    private size_t num_;
    private Mutex mutex_;
    private PileupRunner runner_;

    alias ElementType!(Unqual!(ChunkRange)) Chunk;

    this(string tmp_dir, ChunkRange chunks, BamReader bam, PileupRunner runner) 
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
                        BamReader bam,
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
                                 BamReader bam, PileupRunner runner) {
    return new ChunkDispatcher!ChunkRange(tmp_dir, chunks, bam, runner);
}

int main(string[] args) {
    immutable n_threads = 8;
    auto task_pool = new TaskPool(n_threads);
    auto bam = new BamReader(args[1], task_pool);
    auto ref_file = args[2];
    auto sz = args.length > 3 ? to!size_t(args[3]) : 4_000_000;
    scope(exit) task_pool.finish();

    char[] buf = "/tmp/sambamba-fork-XXXXXX\0".dup;
    mkdtemp(buf.ptr);

    string tmp_dir = to!string(buf.ptr);
//    scope(exit) rmdirRecurse(tmp_dir);

    auto pileupRunner = new PileupRunner(task_pool, std.stdio.stdout, ref_file);

    auto chunks = bam.reads().pileupChunks(false, sz);
    auto dispatcher = chunkDispatcher(tmp_dir, chunks, bam, pileupRunner);
    
    auto threads = new ThreadGroup();
    foreach (i; 0 .. n_threads)
        threads.create(() { worker(dispatcher, bam, task_pool); });

    threads.joinAll();
    pileupRunner.wait();

    return 0;
}
