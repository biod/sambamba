import bamfile;

import std.datetime;
import std.stdio;

void main(string[] args) {
    StopWatch sw;
    sw.start();
    auto bam = BamFile(args[1]);
    foreach (alignment; bam.alignments) {
//        foreach (k, v; alignment.tags) {
//        }
    }
    sw.stop();
    writeln("total time: ", sw.peek().nsecs, "ns");

}
