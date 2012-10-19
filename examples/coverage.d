import bamfile;
import pileuprange;

import std.algorithm, std.c.stdio, std.range;

void main(string[] args) {

    auto dg = delegate (lazy float percentage) {
        static int num;
        num += 1;
        if (num % 8192 == 0) {
            fprintf(stderr, "\r%d%% completed", cast(int)(percentage * 100.0));
        }
    };

    auto bam = BamFile(args[1]);

    auto ref_id = bam.alignments.front.ref_id;
    foreach (column; filter!"a.coverage >= 10"(
                        makePileup(bam.alignmentsWithProgress(dg))))
    {
        printf("%d %d\n", cast(int)column.position, cast(int)column.coverage);
    }

    fprintf(stderr, "\r100%% completed");
}
