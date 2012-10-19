import bamfile;
import pileuprange;

import std.stdio;
import std.datetime;

void main(string[] args) {

    auto bam = BamFile(args[1]);

    StopWatch sw;
    sw.start();
    foreach (read; bam.alignments) {
    }
    sw.stop();
    writeln("Looping through reads: ", sw.peek().msecs, "ms");
    sw.reset();

    sw.start();
    foreach (column; pileupColumns(bam.alignments, false)) {
    }
    sw.stop();
    writeln("Looping through columns (without reference bases): ", sw.peek().msecs, "ms");
    sw.reset();

    sw.start();
    foreach (column; pileupColumns(bam.alignments, true)) {
    }
    sw.stop();
    writeln("Looping through columns (with reference bases): ", sw.peek().msecs, "ms");
}
