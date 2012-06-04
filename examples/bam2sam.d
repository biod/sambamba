import bamfile;
import sam.serialize;

import std.c.stdio;

void main(string[] args) {

    /* This is something like 'samtools view'  */

    auto bam = BamFile(args[1]);

    foreach (alignment; bam.alignments) {
        serialize(alignment, bam.reference_sequences, stdout);
        putchar('\n');
    }
}
