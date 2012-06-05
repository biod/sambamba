import bamfile;

import sam.serialize;
import utils.format;

import std.conv;
import std.c.stdio;
import std.string;
import std.array;

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
        printf(toStringz("usage: " ~ args[0] ~ " <input.bam> <chromosome> <begin> <end>\n"));
        return;
    }

    auto buffer = appender!(ubyte[])();
    buffer.reserve(65536);

    foreach (alignment; bam[chr][beg .. end]) {
        serialize(alignment, bam.reference_sequences, buffer);
        putstring(stdout, cast(string)buffer.data);
        putchar('\n');

        buffer.clear();
    }
}
