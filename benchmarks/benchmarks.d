import bgzfrange;
import bamfile;
import alignment;
import chunkinputstream;

import std.stream;
import std.algorithm;
import std.system;

version(parallel) {
    import std.parallelism;
    TaskPool task_pool;
}

auto bgzfRange(string filename) {
    auto file = new BufferedFile(filename);
    auto compressed_stream = new EndianStream(file, Endian.littleEndian);
    return new BgzfRange(compressed_stream);
}

auto decompressedRange(string filename) {
    version(parallel) {
        return task_pool.map!decompress(bgzfRange(filename), 64);
    } else {
        return map!decompress(bgzfRange(filename));
    }
}

auto unparsedAlignmentRange(string filename) {

    Stream decompressed_stream = makeChunkInputStream(decompressedRange(filename));
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

    return unparsedAlignments(bam);
}

auto parsedAlignmentRange(string filename) {
    version(parallel) {
        return map!parseAlignment(unparsedAlignmentRange(filename));
//        return task_pool.map!parseAlignment(unparsedAlignmentRange(filename), 81920, 10240);
    } else {
        return map!parseAlignment(unparsedAlignmentRange(filename));
    }
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
    version(parallel) {
        auto n_threads = args.length > 1 ? to!uint(args[1]) : totalCPUs;
        task_pool = new TaskPool(n_threads);
        scope(exit) task_pool.finish();
        write(n_threads);
    } else {
        write(1);
    }
    string filename = "../../HG00476.chrom11.ILLUMINA.bwa.CHS.low_coverage.20111114.bam";
    measure!("iterating BGZF blocks", bgzfRange)(filename);
    measure!("iterating decompressed BGZF blocks", decompressedRange)(filename);
    measure!("iterating unparsed alignments", unparsedAlignmentRange)(filename);
    measure!("iterating parsed alignments", parsedAlignmentRange)(filename);
    write("\n");
}
