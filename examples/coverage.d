import bamfile;
import pileuprange;

import std.algorithm, std.c.stdio, std.range;

void main(string[] args) {
    auto bam = BamFile(args[1]);

    auto ref_id = bam.alignments.front.ref_id;
    foreach (column; filter!"a.coverage >= 300"(
                        pileup(
                            filter!"!a.is_unmapped"(
                                until!((Alignment read) { return read.ref_id != ref_id; })(
                                    bam.alignments)))))
    {
        if (column.coverage > 0) {
            printf("%d %d\n", cast(int)column.position, cast(int)column.coverage);
        }
    }
}
