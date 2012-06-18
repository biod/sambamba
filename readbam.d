import bamfile;

import jsonserialization;

import std.datetime;
import std.stdio;

void main(string[] args) {
    StopWatch sw;
    sw.start();
    auto bam = BamFile(args[1]);
    auto count = 0;

    FILE* fp = cast(FILE*)stdout.getFP();
    foreach (alignment; bam.alignments) {
        jsonSerialize(alignment, bam.reference_sequences, fp);
        fputc('\n', fp);
    }
    sw.stop();
}
