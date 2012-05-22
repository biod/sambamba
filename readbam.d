import bamfile;

import std.datetime;
import std.stdio;

void main(string[] args) {
    StopWatch sw;
    sw.start();
    import std.parallelism;
    import std.conv;
    auto tp = new TaskPool(to!int(args[2]) - 1);
    scope(exit) tp.finish();
    auto bam = BamFile(args[1], tp);
    foreach (alignment; bam.alignments) {
//        foreach (k, v; alignment.tags) {
//        }
    }
    sw.stop();
    writeln("total time: ", sw.peek().nsecs, "ns");

}
