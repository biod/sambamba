#!/usr/bin/env rdmd
import bamfile;
import std.stdio;
import std.conv;
import std.range;
  
void main(string[] args) {
    auto bam = BamFile(args[1]);
    if (args.length > 2)
        bam.setBufferSize(to!int(args[2]));
    auto alignments = bam.alignments;
    writeln(walkLength(alignments));
}
