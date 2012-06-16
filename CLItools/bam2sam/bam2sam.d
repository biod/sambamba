import bamfile;
import sam.serialize;
import utils.format;

import std.stdio;
import std.c.stdio : stdout;

void main(string[] args) {

    if (args.length != 2) {
        writeln("usage: " ~ args[0] ~ " <input.bam>");
        writeln("\tprints alignments in SAM format to standard output");
        return;
    }

    BamFile bam;
    try {
        bam = BamFile(args[1]);
    } catch (Throwable e) {
        writeln(e.msg);
        return;
    }

    foreach (alignment; bam.alignments) {
        serialize(alignment, bam.reference_sequences, stdout);
        putcharacter(stdout, '\n');
    }
}
