import bamfile;

import std.datetime;
import std.stdio;

void main(string[] args) {
    auto bam = BamFile(args[1]);

    int[string] tags;
    char[string] typeids;

    foreach (alignment; bam.alignments) {
        foreach (key, value; alignment) {
            tags[key] += 1;
            typeids[key] = value.bam_typeid;
        }
    }

    foreach (key, count; tags) {
        writeln("[", key, "] - ", count, " (", typeids[key], ")");
    }
}
