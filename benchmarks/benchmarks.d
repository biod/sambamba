import bgzfrange;
import bamfile;
import alignment;
import alignmentrange;
import chunkinputstream;

import std.stream;
import std.algorithm;
import std.system;

version(parallel) {
    import std.parallelism;
    int n_threads;
    TaskPool task_pool;
}

auto bgzfRange(string filename) {
    auto file = new BufferedFile(filename);
    auto compressed_stream = new EndianStream(file, Endian.littleEndian);
    return BgzfRange(compressed_stream);
}

auto decompressedRange(string filename) {
    version(parallel) {
        return task_pool.map!decompressBgzfBlock(bgzfRange(filename), 64);
    } else {
        return map!decompressBgzfBlock(bgzfRange(filename));
    }
}

auto getAlignmentRange(string filename) {

    IChunkInputStream decompressed_stream = makeChunkInputStream(decompressedRange(filename));
    Stream bam = new EndianStream(decompressed_stream, Endian.littleEndian); 

    bam.readString(4); // skip magic
    int l_text;
    bam.read(l_text);
    bam.readString(l_text); // skip header
    int n_ref;
    bam.read(n_ref);
    while (n_ref-- > 0) {
        int l_name;
        bam.read(l_name);
        bam.readString(l_name);
        int l_ref;
        bam.read(l_ref);
    } // skip reference sequences information

    return alignmentRange(decompressed_stream);
}

import std.datetime;
import std.stdio;
import std.conv;

ulong measure(string desc, alias func, Args...)(Args args) {
    StopWatch sw;
    sw.start();
    auto range = func(args);
    foreach (elem; range) {
    }
    sw.stop();
    write("\t", sw.peek().msecs);
    return sw.peek().nsecs;
}

void main(string[] args) {
    if (args.length == 1) {
        writeln("usage: " ~ args[0] ~ " <input.bam> <number of threads>");
        writeln("    prints the following numbers:");
        writeln("    - number of threads used for decompression");
        writeln("    - time to iterate over BGZF blocks");
        writeln("    - time to iterate over decompressed blocks");
        writeln("    - time to iterate over alignments");
        writeln();
        writeln("    For serial version, 1 thread is used.");
        writeln();
        writeln("    Also, timings for first step depend a lot on");
        writeln("    whether the file is cached in memory or not.");
        return;
    }
    version(parallel) {
        n_threads = args.length > 2 ? to!uint(args[2]) : totalCPUs;
        task_pool = new TaskPool(n_threads);
        scope(exit) task_pool.finish();
        write(n_threads);
    } else {
        write(1);
    }
    string filename = args[1];
    measure!("iterating BGZF blocks", bgzfRange)(filename);
    measure!("iterating decompressed BGZF blocks", decompressedRange)(filename);
    measure!("iterating alignments", getAlignmentRange)(filename);
    write("\n");
}
