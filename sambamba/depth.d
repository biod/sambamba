/*
    This file is part of Sambamba.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.depth;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.region;
import bio.bam.multireader;
import sambamba.utils.common.bed;
import sambamba.utils.common.filtering;
import sambamba.utils.common.intervaltree;

import std.stdio;
import std.exception;
import std.ascii: isWhite;
import std.range;
import std.algorithm;
import std.array;
import std.getopt;
import std.parallelism;
import std.typecons;

version(standalone) {
    int main(string[] args) {
        return depth_main(args);
    }
}

alias uint pos_t;

void printUsage() {
    stderr.writeln("Usage: sambamba-depth [options] input.bam  [input2.bam [...]]");
    stderr.writeln();
    stderr.writeln("          All BAM files must be coordinate-sorted and indexed.");
    stderr.writeln();
    stderr.writeln("Options: -F, --filter=FILTER");
    stderr.writeln("                    set custom filter for alignments");
    stderr.writeln("         -L, --regions=FILENAME");
    stderr.writeln("                    list or regions of interest");  
    stderr.writeln("         -b, --per-base-output-fn");
    stderr.writeln("                    specify output filename for per-base output");
    stderr.writeln("         -r, --per-region-output-fn");
    stderr.writeln("                    specify output filename for per-region output");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
}

struct RegionStats {
    size_t id;
    private bool first_occurrence = true;

    size_t n_bases; // all bases
    size_t n_good_bases; // determined by  mapping quality / base quality threshold

    size_t n_reads; // all reads overlapping
    size_t n_good_reads;

    size_t percentage_covered; // 
    size_t percentage_good_covered; // 
}

abstract class RegionStatsCollector {
    void nextColumn(uint ref_id, pos_t position,
                    scope void delegate(ref RegionStats) updater);

    const(RegionStats)[] regionStatistics();
}

class GeneralRegionStatsCollector : RegionStatsCollector {
    private BamRegion[] bed_;

    private {
        // stores region statistics and index of the region in initial array
        alias Tuple!(RegionStats, size_t) payload;

        alias IntervalTree!(payload, uint) intervalTree;
        alias IntervalTreeNode!(payload, uint) intervalTreeNode;
        intervalTree[] trees_;
    }

    this(BamRegion[] bed) {
        bed_ = bed;

        if (bed.length == 0)
            return;

        trees_.length = bed_.back.ref_id + 1;

        size_t start_index = 0;
        size_t end_index = start_index;

        intervalTreeNode[] intervals;

        while (start_index < bed.length) {
            while (end_index < bed.length && bed[end_index].ref_id == bed[start_index].ref_id)
                ++end_index;

            intervals.length = end_index - start_index;
            foreach (i; 0 .. intervals.length) {
                auto start = bed[start_index + i].start;
                auto stop = bed[start_index + i].end;
                RegionStats stats;
                stats.id = start_index + i;
                auto value = tuple(stats, start_index + i);
                intervals[i] = new intervalTreeNode(start, stop, value);
            }

            trees_[bed[start_index].ref_id] = new intervalTree(intervals);
            start_index = end_index;
        }
    }

    override void nextColumn(uint ref_id, pos_t position,
                             scope void delegate(ref RegionStats) updater)
    {
        if (ref_id >= trees_.length || trees_[ref_id] is null)
            return;

        foreach (node; trees_[ref_id].eachOverlap(position, position + 1)) {
            updater(node.value[0]);
        }
    }

    override const(RegionStats)[] regionStatistics() {
        RegionStats[] reg_stats;
        reg_stats.length = bed_.length;
        foreach (i; 0 .. trees_.length) {
            if (trees_[i] !is null) {
                foreach (node; trees_[i].eachOverlap(pos_t.min, pos_t.max))
                    reg_stats[node.value[1]] = node.value[0];
            }
        }
        return reg_stats;
    }
}

bool isSortedAndNonOverlapping(BamRegion[] bed) {
    size_t n = bed.length;
    if (n <= 1) return true;
    foreach (k; 0 .. n - 1) {
        auto reg1 = bed[k];
        auto reg2 = bed[k + 1];
        if (reg1.ref_id > reg2.ref_id)
            return false;
        if (reg1.ref_id < reg2.ref_id)
            continue;
        if (reg1.end > reg2.start)
            return false;
    }
    return true;
}

class NonOverlappingRegionStatsCollector : RegionStatsCollector {
    private BamRegion[] bed_;
    private RegionStats[] region_statistics_;

    size_t current_index_;

    this(BamRegion[] bed) {
        assert(isSortedAndNonOverlapping(bed));
        bed_ = bed;
        region_statistics_.length = bed_.length;
        current_index_ = 0;
    }

    override void nextColumn(uint ref_id, pos_t position,
                             scope void delegate(ref RegionStats) updater)
    {
        while (bed_.length > 0 &&
               bed_.front.fullyLeftOf(ref_id, position))
        {
            bed_ = bed_[1 .. $];
            ++current_index_;
        }

        if (bed_.length > 0 &&
            bed_.front.overlaps(ref_id, position.to!uint))
        {
            updater(region_statistics_[current_index_]);
        }
    }

    override const(RegionStats)[] regionStatistics() {
        return region_statistics_;
    }
}

int depth_main(string[] args) {

    string bed_filename = null;
    string base_fn = null;
    string region_fn = null;
    string query = null;

    File per_base_output = stdout;
    File per_region_output;

    BamRegion[] raw_bed; // may be overlapping
    string[] raw_bed_lines;

    int n_threads;
    
    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               "filter|F",               &query,
               "regions|L",              &bed_filename,
               "per-base-output-fn|b",   &base_fn,
               "per-region-output-fn|r", &region_fn,
               "nthreads|t",             &n_threads);

        if (args.length < 2) {
            printUsage();
            return 0;
        }

        defaultPoolThreads = max(n_threads, 0);

        Filter filter = new NullFilter();
        if (query !is null) {
            filter = createFilterFromQuery(query);
        }

        auto bam_filenames = args[1 .. $];
        auto bam = new MultiBamReader(bam_filenames);
        enforce(bam.header.sorting_order == SortingOrder.coordinate,
                "All files must be coordinate-sorted");
        enforce(bam.has_index, "All files must be indexed");

        InputRange!(MultiBamRead!BamRead) reads;
        if (bed_filename !is null) {
            auto bed = parseBed(bed_filename, bam);
            raw_bed = parseBed(bed_filename, bam, false, &raw_bed_lines);
            reads = inputRangeObject(bam.getReadsOverlapping(bed));
        } else {
            reads = inputRangeObject(bam.reads);
        }

        auto filtered_reads = filtered(reads, filter);
        auto pileup = makePileup(filtered_reads);
        auto cov = new uint[bam_filenames.length];

        size_t current_region_index = 0;

        size_t[] overlapping_region_indices;

        RegionStatsCollector stats_collector = null;
        if (bed_filename !is null && region_fn !is null) {
            if (isSortedAndNonOverlapping(raw_bed))
                stats_collector = new NonOverlappingRegionStatsCollector(raw_bed);
            else
                stats_collector = new GeneralRegionStatsCollector(raw_bed);

            per_region_output = File(region_fn, "w+");
        }

        if (base_fn) {
            per_base_output = File(base_fn, "w+");
        }

        foreach (column; pileup) {
            auto ref_name = bam.reference_sequences[column.ref_id].name;
            foreach (read; column.reads)
                cov[read.file_id] += 1;

            with (per_base_output) {
                write(ref_name, '\t', column.position);
                foreach (partial_cov; cov)
                    write("\t", partial_cov);
                writeln();
            }

            cov[] = 0;

            if (stats_collector !is null && column.ref_id >= 0) {
                stats_collector.nextColumn(column.ref_id.to!uint,
                                           column.position.to!pos_t,
                (ref RegionStats stats) {
                    with (stats) {
                        n_bases += column.coverage;

                        if (first_occurrence) {
                            n_reads += column.coverage;
                            first_occurrence = false;
                        } else {
                            n_reads += column.reads_starting_here.length;
                        }
                    }});
            }
        }

        if (stats_collector !is null) {
            auto region_statistics = stats_collector.regionStatistics();
            assert(region_statistics.length == raw_bed.length);
            foreach (i, region; raw_bed) {
                auto stats = region_statistics[i];
                with(per_region_output) {
                    write(raw_bed_lines[i]);
                    // writeln('\t', stats.n_bases, '\t', stats.n_reads);
                    if (!isWhite(raw_bed_lines[i].back))
                        write('\t');
                    writeln(stats.n_reads);
                }
            }
        }

        return 0;

    } catch (Exception e) {
        stderr.writeln("sambamba-depth: ", e.msg);

        version(development) {
            throw e;
        }
        return 1;
    }

    return 0;
}
