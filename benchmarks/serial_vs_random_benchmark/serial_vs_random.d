import bamfile;

import std.stdio;
import std.datetime;
import std.range;

void printUsage() {
    stderr.writeln(
"usage: ./iterate <input.bam> METHOD

     Utility for measuring time that it takes to iterate the whole file using METHOD.
     All reads in <input.bam> must be aligned to the same reference sequence.

     Prints number of alignments and time it took to iterate them.

     METHOD can be either SERIAL or RANDOM. The meaning is as follows:

        SERIAL is when the file is read serially from beginning to end
               (that's what happens when BamFile.alignments is used)

        RANDOM is when degenerate case of random access is used, i.e.
               when region ref:0-length(ref) is retrieved using the same
               procedure as used for fetching arbitrary region. That is,
               this includes associated overhead of using index file,
               maintaining a cache (which, in this particular case, is not
               needed at all).


     The idea is that even degenerate random access is quite close in performance to
     serial reading, despite a significant amount of introduced overhead.
     ");
}

void main(string[] args) {
    if (args.length < 3) {
        printUsage();
        return;
    }

    auto bam = BamFile(args[1]);

    switch (args[2]) {
        case "SERIAL":
            StopWatch sw;
            sw.start();
            writeln(walkLength(bam.alignments));
            sw.stop();
            writeln(sw.peek().msecs, "ms");
            return;
        case "RANDOM":
            auto ref_id = bam.alignments.front.ref_id;
            StopWatch sw;
            sw.start();
            writeln(walkLength(bam.reference(ref_id)[]));
            sw.stop();
            writeln(sw.peek().msecs, "ms");
            return;
        default:
            printUsage();
            return;
    }
}
