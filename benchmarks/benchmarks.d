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
}

auto bgzfRange(string filename) {
    auto file = new BufferedFile(filename);
    auto compressed_stream = new EndianStream(file, Endian.littleEndian);
    return BgzfRange(compressed_stream);
}

auto decompressedRange(string filename) {
    version(parallel) {
//        return task_pool.map!decompressBgzfBlock(bgzfRange(filename), 64);
        import utils.range;
        return parallelTransform!decompressBgzfBlock(bgzfRange(filename), 64, n_threads);
    } else {
        return map!decompressBgzfBlock(bgzfRange(filename));
    }
}

auto unparsedAlignmentRange(string filename) {

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

    return unparsedAlignments!(IteratePolicy.withoutOffsets)(decompressed_stream);
}

auto parsedAlignmentRange(string filename) {
    return map!makeAlignment(unparsedAlignmentRange(filename));
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
        n_threads = args.length > 1 ? to!uint(args[1]) : totalCPUs;
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
