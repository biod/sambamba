import bamfile;

import std.datetime;
import std.stdio;

void main(string[] args) {
    StopWatch sw;
    sw.start();
    auto bam = BamFile(args[1]);
    auto count = 0;
    foreach (alignment; bam.alignments) {
        count += 1;
    }
    sw.stop();
    writeln("total time: ", sw.peek().nsecs, "ns");
    writeln("alignments found: ", count);

}
