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
module sambamba.view;

import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.region;
import bio.bam.multireader;
import sambamba.utils.common.bed;
import sambamba.utils.common.filtering;

import std.stdio;
import std.exception;
import std.range;
import std.algorithm;
import std.array;
import std.getopt;
import std.parallelism;

version(standalone) {
    int main(string[] args) {
        return depth_main(args);
    }
}

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
    size_t n_bases; // all bases
    size_t n_good_bases; // determined by  mapping quality / base quality threshold

    size_t percentage_covered; // 
    size_t percentage_good_covered; // 
}

int depth_main(string[] args) {

    string bed_filename;
    string base_fn;
    string region_fn;
    string query;

    BedIndex intervals; // may include overlapping intervals

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
        if (bed_filename !is null) {
            intervals = readIntervals(bed_filename, false);
        }

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
            reads = inputRangeObject(bam.getReadsOverlapping(bed));
        } else {
            reads = inputRangeObject(bam.reads);
        }

        auto filtered_reads = filtered(reads, filter);
        auto pileup = makePileup(filtered_reads);
        auto cov = new uint[bam_filenames.length];
        foreach (column; pileup) {
            auto ref_name = bam.reference_sequences[column.ref_id].name;
            foreach (read; column.reads)
                cov[read.file_id] += 1;
            write(ref_name, '\t', column.position);
            foreach (partial_cov; cov)
                write("\t", partial_cov);
            writeln();
            cov[] = 0;
        }
        return 0;

    } catch (Exception e) {
        stderr.writeln("sambamba-depth: ", e.msg);

        version(development) {
            throw e;
        }
        return 1;
    }
}
