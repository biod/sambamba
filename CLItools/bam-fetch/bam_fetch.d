import bamfile;

import sam.serialize;

import std.conv;
import std.c.stdio;
import std.string;

void main(string[] args) {
    BamFile bam;
    string chr;
    int beg;
    int end;
    try {
        bam = BamFile(args[1]);
        chr = args[2];
        beg = to!int(args[3]);
        end = to!int(args[4]);
    } catch (Throwable e) {
        printf(toStringz("usage: " ~ args[0] ~ " <input.bam> <chromosome> <begin> <end>"));
        return;
    }

    foreach (alignment; bam[chr][beg .. end]) {
        serialize(alignment, bam.reference_sequences, stdout);
        putchar('\n');
    }
}
