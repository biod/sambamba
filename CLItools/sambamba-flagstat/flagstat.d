/// port of samtools flagstat tool
import bamfile;
import std.stdio, std.conv, std.parallelism;

ulong[2] reads, pair_all, pair_good, first, second, single, pair_map, mapped,
         dup, diff_chr, diff_high;

void writeParam(string description, ulong[2] param) {
    writefln("%s + %s %s", param[0], param[1], description);
}

float percent(ulong a, ulong b) { return to!float(a) / b * 100.0; }

void writeParamWithPercentage(string description, ulong[2] param, ulong[2] total) {
    writefln("%s + %s %s (%.2f%%:%.2f%%)", param[0], param[1], description,
             percent(param[0], total[0]), percent(param[1], total[1]));
}

int main(string[] args) {
    if (args.length == 1 || args.length > 3) {
        stderr.writeln("Usage: sambamba-flagstat <input.bam> [nthreads=#cores]");
        return 1;
    }

    try {
    auto threads = args.length == 2 ? totalCPUs : to!uint(args[2]);
    auto task_pool = new TaskPool(threads);
    scope(exit) task_pool.finish();

    auto bam = BamFile(args[1], task_pool);

    foreach (read; bam.alignments) {
        size_t failed = read.failed_quality_control ? 1 : 0;
        ++reads[failed];
        if (!read.is_unmapped) ++mapped[failed];
        if (read.is_duplicate) ++dup[failed];
        if (read.is_paired) {
            ++pair_all[failed];
            if (read.proper_pair) ++pair_good[failed];
            if (read.is_first_of_pair) ++first[failed];
            if (read.is_second_of_pair) ++second[failed];
            if (read.mate_is_unmapped && !read.is_unmapped) ++single[failed];
            if (!read.is_unmapped && !read.mate_is_unmapped) {
                ++pair_map[failed];
                if (read.ref_id != read.next_ref_id) {
                    ++diff_chr[failed];
                    if (read.mapping_quality >= 5)
                        ++diff_high[failed];
                }
            }
        }
    }

    scope(exit) {
        writeParam("in total (QC-passed reads + QC-failed reads)", reads);
        writeParam("duplicates", dup);
        writeParamWithPercentage("mapped", mapped, reads);
        writeParam("paired in sequencing", pair_all);
        writeParam("read1", first);
        writeParam("read2", second);
        writeParamWithPercentage("properly paired", pair_good, pair_all);
        writeParam("with itself and mate mapped", pair_map);
        writeParamWithPercentage("singletons", single, pair_all);
        writeParam("with mate mapped to a different chr", diff_chr);
        writeParam("with mate mapped to a different chr (mapQ>=5)", diff_high);
    }
    } catch (Throwable e) {
        stderr.writeln(e.msg);
    }
    return 0;
}
