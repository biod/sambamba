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
import bio.core.base;
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
    stderr.writeln("Usage: sambamba-depth region|window [options] input.bam  [input2.bam [...]]");
    stderr.writeln();
    stderr.writeln("          All BAM files must be coordinate-sorted and indexed.");
    stderr.writeln();
    stderr.writeln("          The tool has three modes: base, region, and window,");
    stderr.writeln("          each name means per which unit to print the statistics.");
    stderr.writeln();
    stderr.writeln("Common options:");
    stderr.writeln("         -F, --filter=FILTER");
    stderr.writeln("                    set custom filter for alignments");
    stderr.writeln("         -o, --output-prefix=PREFIX");
    stderr.writeln("                    output filename prefix");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -c, --min-coverage=MAXCOVERAGE");
    stderr.writeln("                    minimum mean coverage for output");
    stderr.writeln("         -C, --max-coverage=MAXCOVERAGE");
    stderr.writeln("                    maximum mean coverage for output");
    stderr.writeln("         -a, --annotate");
    stderr.writeln("                    add additional column of y/n instead of");
    stderr.writeln("                    skipping records not satisfying the criteria");
    stderr.writeln("base subcommand options:");
    stderr.writeln("         -L, --regions=FILENAME");
    stderr.writeln("                    list or regions of interest (optional)");
    stderr.writeln("         -q, --min-base-quality=QUAL");
    stderr.writeln("                    don't count bases with lower base quality");
    stderr.writeln("region subcommand options:");
    stderr.writeln("         -L, --regions=FILENAME");
    stderr.writeln("                    list or regions of interest (required)");
    stderr.writeln("         --cov-threshold=COVTHRESHOLD");
    stderr.writeln("                    multiple thresholds can be provided,");
    stderr.writeln("                    for each one an extra column will be added,");
    stderr.writeln("                    the percentage of bases in the region");
    stderr.writeln("                    where coverage is more than this value");
    stderr.writeln("         --combined");
    stderr.writeln("                    output single file with combined statistics");
    stderr.writeln("                    including all samples");
    stderr.writeln("window subcommand options:");
    stderr.writeln("         --window-size=WINDOWSIZE");
    stderr.writeln("                    breadth of the window, in bp (required)");
    stderr.writeln("         --overlap=OVERLAP");
    stderr.writeln("                    overlap of successive windows, in bp (default is 0)");
    stderr.writeln("         --cov-threshold=COVTHRESHOLD,");
    stderr.writeln("         --combined");
    stderr.writeln("                    same meaning as in 'region' subcommand");
}

enum Mode {
    base,
    region,
    window
}

abstract class RegionStatsCollector {
    void nextColumn(uint ref_id, pos_t position,
                    scope void delegate(size_t id) updater);
}

class GeneralRegionStatsCollector : RegionStatsCollector {
    private BamRegion[] bed_;

    private {
        // stores index of the region in initial array
        alias size_t payload;
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
                auto value = start_index + i;
                intervals[i] = new intervalTreeNode(start, stop, value);
            }

            trees_[bed[start_index].ref_id] = new intervalTree(intervals);
            start_index = end_index;
        }
    }

    override void nextColumn(uint ref_id, pos_t position,
                             scope void delegate(size_t id) updater)
    {
        if (ref_id >= trees_.length || trees_[ref_id] is null)
            return;

        foreach (node; trees_[ref_id].eachOverlap(position, position + 1)) {
            updater(node.value);
        }
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

    size_t current_index_;

    this(BamRegion[] bed) {
        assert(isSortedAndNonOverlapping(bed));
        bed_ = bed;
        current_index_ = 0;
    }

    override void nextColumn(uint ref_id, pos_t position,
                             scope void delegate(size_t id) updater)
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
            updater(current_index_);
        }
    }
}

class WindowStatsCollector : RegionStatsCollector {
    MultiBamReader bam;
    size_t window_size;
    size_t overlap;
    size_t step;
    size_t n;

    this(MultiBamReader bam, size_t window_size, size_t overlap, size_t n) {
        this.bam = bam;
        this.window_size = window_size;
        this.overlap = overlap;
        step = window_size - overlap;
        this.n = n;
    }

    override void nextColumn(uint ref_id, pos_t position,
                             scope void delegate(size_t id) updater)
    {
        size_t k;
        if (position < window_size)
            k = position / step + 1;
        else
            k = n;

        foreach (id; 0 .. k)
            updater(id);
    }
}

struct CustomBamRead {
    MultiBamRead!BamRead read;
    alias read this;

    this(MultiBamRead!BamRead read, string[string] rg2sm) {
        this.read = read;
        auto rg_value = read["RG"];
        if (!rg_value.is_nothing) {
            auto rg_str = *(cast(string*)(&rg_value));
            sample_name = rg2sm[rg_str];
        }
    }

    string sample_name;

    CustomBamRead dup() @property const {
        CustomBamRead r = void;
        r.read = read.dup;
        r.sample_name = sample_name;
        return r;
    }
}

alias Column = PileupRange!(InputRange!CustomBamRead).Column;

abstract class ColumnPrinter {
    double min_cov = 0.0;
    double max_cov = 1e50;
    bool annotate = false;

    string output_prefix;
    MultiBamReader bam;

    BamRegion[] raw_bed;

    void setBed(BamRegion[] bed) {
        raw_bed = bed;
    }

    string[] raw_bed_lines;

    abstract void init(ref string[] args);
    abstract void push(ref Column);
    abstract void close();
}

class PerBasePrinter : ColumnPrinter {
    ubyte min_base_quality;
    NonOverlappingRegionStatsCollector stats_collector;
    private File _f;

    override void init(ref string[] args) {
        getopt(args,
               std.getopt.config.caseSensitive,
               "min-base-quality|q", &min_base_quality);

	if (output_prefix is null)
	    _f = stdout;
	else 
	    _f = File(output_prefix ~ ".cov");

	_f.writeln("REF\tPOS\tCOV\tA\tC\tG\tT\tDEL\tREFSKIP");
    }

    override void setBed(BamRegion[] bed) {
	super.setBed(bed);
	stats_collector = new NonOverlappingRegionStatsCollector(bed);
    }

    private void writeColumn(ref File f, ref Column c) {
	size_t[5] coverage;
	size_t deletions;
	size_t ref_skips;
	foreach (read; c.reads) {
	    if (read.current_base == '-') {
		if (read.cigar_operation.type == 'D')
		    deletions += 1;
		else
		    ref_skips += 1;
		continue;
	    }
		
	    if (read.current_base_quality >= min_base_quality)
		coverage[Base5(read.current_base).internal_code] += 1;
	}

	size_t total_coverage = reduce!`a+b`(coverage[]) + deletions + ref_skips;

	bool ok = total_coverage >= min_cov || total_coverage <= max_cov;
	if (!ok && !annotate)
	    return;
	
	_f.write(bam.reference_sequences[c.ref_id].name, '\t',
		 c.position, '\t', total_coverage);
	foreach (i; 0 .. 4)
	    _f.write('\t', coverage[i]);
	_f.write('\t', deletions, '\t', ref_skips);

	if (annotate)
	    _f.write('\t', ok ? 'y' : 'n');
	_f.writeln();
    }

    override void push(ref Column c) {
	if (stats_collector !is null) {
	    bool output = false;
	    stats_collector.nextColumn(cast(uint)c.ref_id, cast(uint)c.position,
				       (size_t id) { output = true; });
	    if (output)
		writeColumn(_f, c); 
	} else {
	    writeColumn(_f, c);
	}
    }

    override void close() {
    }
}

class PerSampleRegionData {
    this(size_t n_coverage_counters, size_t n_regions) {
        coverage_counters_ = new uint[][](n_coverage_counters, n_regions);
        n_reads_.length = n_regions;
        n_bases_.length = n_regions;
    }

    // for each coverage threshold, we hold here numbers of bases with
    // cov. >= that threshold, for each region
    private uint[][] coverage_counters_;
    private uint[] n_reads_; // for each region
    private uint[] n_bases_; // ditto

    ref uint coverage_count(size_t cov_id, size_t region_id) {
        return coverage_counters_[cov_id][region_id];
    }

    ref uint n_reads(size_t id) { return n_reads_[id]; }
    ref uint n_bases(size_t id) { return n_bases_[id]; }

    void reset(size_t id) {
        n_reads_[id] = 0;
        n_bases_[id] = 0;
        foreach (ref cov_counter; coverage_counters_)
            cov_counter[id] = 0;
    }
}

private size_t overlap(R)(BamRegion region, auto ref R read) {
    assert(region.ref_id == read.ref_id);
    size_t s1 = region.start;
    size_t e1 = region.end;
    size_t s2 = read.position;
    size_t e2 = read.end_position;
    if (e1 <= s2 || e2 <= s1)
        return 0;
    return min(e1, e2) - max(s1, s2);
}

abstract class PerRegionPrinter : ColumnPrinter {
    RegionStatsCollector stats_collector;
    private PerSampleRegionData[string] samples;

    string outputFilename(string sample_name) {
        return output_prefix ~ sample_name ~ ".bed";
    }

    private string getSampleName(R)(auto ref R read) {
        if (combined)
            return "";
        else
            return read.sample_name;
    }

    private void countRead(R)(auto ref R read, size_t id) {
        auto sample = getSampleName(read);
        auto data = getSampleData(sample);
        data.n_reads(id) += 1;
        data.n_bases(id) += overlap(getRegionById(id), read);
    }

    abstract BamRegion getRegionById(size_t id);
    abstract PerSampleRegionData getSampleData(string sample);
    abstract bool isFirstOccurrence(size_t id);
    abstract void markAsSeen(size_t id);
    abstract void writeOriginalBedLine(ref File f, size_t id);

    bool combined = false;

    uint[] cov_thresholds;
    bool[] is_first_occurrence;

    uint[string] cov_per_sample;

    override void init(ref string[] args) {
        getopt(args,
               std.getopt.config.caseSensitive,
               "cov-threshold", &cov_thresholds,
               "combined", &combined);
    }

    override void push(ref Column column) {
        uint ref_id = column.ref_id.to!uint;
        pos_t position = column.position.to!pos_t;

        stats_collector.nextColumn(ref_id, position,
        (size_t id) {
            if (isFirstOccurrence(id)) {
                foreach (read; column.reads)
                    countRead(read, id);
                markAsSeen(id);
            } else {
                foreach (read; column.reads_starting_here)
                    countRead(read, id);
            }

            foreach (sample, ref cov; cov_per_sample)
                cov = 0;

            foreach (read; column.reads) {
                auto sample = getSampleName(read);
                cov_per_sample[sample] += 1;
            }

            foreach (sample, ref cov; cov_per_sample) {
                auto data = getSampleData(sample);
                foreach (i, threshold; cov_thresholds)
                    if (cov >= threshold)
                        data.coverage_count(i, id) += 1;
            }
        });
    }

    void printRegionStats(ref File f, size_t id, PerSampleRegionData data) {
        auto region = getRegionById(id);
        auto length = region.end - region.start;
        with(f) {
            auto mean_cov = data.n_bases(id).to!float / length;

            bool ok = mean_cov >= this.min_cov && mean_cov <= this.max_cov;

            if (!ok && !this.annotate)
                return;

            writeOriginalBedLine(f, id);
            write(data.n_reads(id), '\t', mean_cov);
            foreach (j; 0 .. cov_thresholds.length)
                write('\t', data.coverage_count(j, id).to!float * 100 / length);

            if (annotate)
                write('\t', ok ? 'y' : 'n');

            writeln();
            flush();
        }
    }
}

class PerBedRegionPrinter : PerRegionPrinter {
    bool[] is_first_occurrence;

    override PerSampleRegionData getSampleData(string sample) {
        auto ptr = sample in samples;
        if (ptr)
            return *ptr;
        auto data = new PerSampleRegionData(cov_thresholds.length,
                                            raw_bed.length);
        samples[sample] = data;
        return data;
    }

    override bool isFirstOccurrence(size_t id) {
        return is_first_occurrence[id];
    }

    override void markAsSeen(size_t id) {
        is_first_occurrence[id] = false;
    }

    override void writeOriginalBedLine(ref File f, size_t id) {
        f.write(raw_bed_lines[id]);
        if (!isWhite(raw_bed_lines[id].back))
            f.write('\t');
    }

    override BamRegion getRegionById(size_t id) {
        return raw_bed[id];
    }

    override void setBed(BamRegion[] bed) {
        raw_bed = bed;
        is_first_occurrence = new bool[](raw_bed.length);
        is_first_occurrence[] = true;

        if (isSortedAndNonOverlapping(raw_bed))
            stats_collector = new NonOverlappingRegionStatsCollector(raw_bed);
        else
            stats_collector = new GeneralRegionStatsCollector(raw_bed);
    }

    override void close() {
        foreach (sample_name, data; samples) {
            auto f = File(outputFilename(sample_name), "w+");
            foreach (id; 0 .. raw_bed.length) {
                printRegionStats(f, id, data);
            }
        }
    }
}

class PerWindowPrinter : PerRegionPrinter {
    size_t window_size;
    size_t overlap;
    size_t step;
    size_t n; // number of windows to keep at each moment

    bool[] is_first_occurrence;

    size_t leftmost_window_index = 0;
    size_t leftmost_window_start_pos = 0;
    int window_ref_id = -1;
    size_t ref_length;

    File[string] output_files;
    ref File outputFile(string sample) {
        auto ptr = sample in output_files;
        if (ptr !is null)
            return *ptr;
        output_files[sample] = File(outputFilename(sample), "w+");
        return output_files[sample];
    }

    void printWindowStats(size_t id) {
        foreach (sample, data; samples) {
            printRegionStats(outputFile(sample), leftmost_window_index, data);
        }
    }

    void resetAllWindows() {
        foreach (sample, data; samples) {
            foreach (id; 0 .. n) {
                data.reset(id);
            }
        }
        is_first_occurrence[] = true;
        leftmost_window_index = 0;
        leftmost_window_start_pos = 0;
    }

    void finishLeftMostWindow() {
        printWindowStats(leftmost_window_index);
        foreach (sample, data; samples)
            data.reset(leftmost_window_index);
        is_first_occurrence[leftmost_window_index] = true;

        leftmost_window_index += 1;
        if (leftmost_window_index == n)
            leftmost_window_index = 0;
        leftmost_window_start_pos += step;
    }

    size_t windowStart(size_t id) {
        size_t k_from_leftmost;
        if (id >= leftmost_window_index) {
            k_from_leftmost = id - leftmost_window_index;
        } else {
            k_from_leftmost = n - leftmost_window_index + id;
        }
        return leftmost_window_start_pos + step * k_from_leftmost;
    }

    override PerSampleRegionData getSampleData(string sample) {
        auto ptr = sample in samples;
        if (ptr)
            return *ptr;
        auto data = new PerSampleRegionData(cov_thresholds.length, n);
        samples[sample] = data;
        return data;
    }

    override bool isFirstOccurrence(size_t id) {
        return is_first_occurrence[id];
    }

    override void markAsSeen(size_t id) {
        is_first_occurrence[id] = false;
    }

    override BamRegion getRegionById(size_t id) {
        auto start = windowStart(id);
        auto end = start + window_size;
        return BamRegion(window_ref_id, cast(uint)start, cast(uint)end);
    }

    override void writeOriginalBedLine(ref File f, size_t id) {
        auto region = getRegionById(id);
        f.write(bam.reference_sequences[region.ref_id].name, '\t',
                region.start, '\t',
                region.end, '\t');
    }

    override void init(ref string[] args) {
        getopt(args, std.getopt.config.caseSensitive,
               "window-size", &window_size,
               "overlap", &overlap,
               "cov-threshold", &cov_thresholds,
               "combined", &combined);

        enforce(window_size > 0,
                "positive window size must be specified");

        enforce(overlap < window_size,
                "specified overlap is larger than window size");

        step = window_size - overlap;
        n = window_size / step;
        if (window_size % step != 0)
            ++n;

        is_first_occurrence = new bool[n];
        is_first_occurrence[] = false;

        stats_collector = new WindowStatsCollector(bam, window_size, overlap, n);
    }

    private void printEmptyWindows(int ref_id) {
        window_ref_id = ref_id;
        foreach (j; 0 .. bam.reference_sequences[ref_id].length / step)
            finishLeftMostWindow();
        resetAllWindows();
    }

    private void moveToReference(int ref_id) {
        window_ref_id = ref_id;
        ref_length = bam.reference_sequences[ref_id].length;
    }

    override void push(ref Column column) {
        if (window_ref_id == -1) {
            foreach (int k; 0 .. column.ref_id)
                printEmptyWindows(k);
            moveToReference(column.ref_id);
        } else if (column.ref_id != window_ref_id) {
            while (leftmost_window_start_pos + window_size < ref_length)
                finishLeftMostWindow();
            resetAllWindows();
            foreach (k; window_ref_id + 1 .. column.ref_id)
                printEmptyWindows(k);
            moveToReference(column.ref_id);
        }

        while (column.position >= leftmost_window_start_pos + window_size)
            finishLeftMostWindow();
        super.push(column);
    }

    override void close() {
        while (leftmost_window_start_pos + window_size < ref_length)
            finishLeftMostWindow();
        foreach (k; window_ref_id + 1 .. bam.reference_sequences.length)
            printEmptyWindows(cast(int)k);
        foreach (sample_name, f; output_files)
            f.close();
    }
}

int depth_main(string[] args) {

    int n_threads;
    string query = null;
    Filter read_filter = new NullFilter();

    string bed_filename = null;

    if (args.length < 3) {
        printUsage();
        return 0;
    }

    Mode mode;
    ColumnPrinter printer;

    switch (args[1]) {
    case "base":
        mode = Mode.base;
        printer = new PerBasePrinter();
        break;
    case "region":
        mode = Mode.region;
        printer = new PerBedRegionPrinter();
        break;
    case "window":
        mode = Mode.window;
        printer = new PerWindowPrinter();
        break;
    default:
        printUsage();
        return 0;
    }

    args = args[1 .. $];

    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               std.getopt.config.passThrough,
               "filter|F",               &query,
               "output-prefix|o",        &printer.output_prefix,
               "nthreads|t",             &n_threads,
               "min-coverage|c",         &printer.min_cov,
               "max-coverage|C",         &printer.max_cov,
               "annotate|a",             &printer.annotate);

        if (mode != Mode.window) {
            getopt(args,
                   std.getopt.config.caseSensitive,
                   std.getopt.config.passThrough,
                   "regions|L", &bed_filename);
        }

	if (mode != Mode.base) {
	    enforce(printer.output_prefix !is null,
		    "Output prefix must be provided for the " ~ mode.to!string ~
		    " mode");
	}

        // handles subcommand arguments and removes them from the list
        printer.init(args);

        defaultPoolThreads = max(n_threads, 0);

        if (query !is null) {
            read_filter = createFilterFromQuery(query);
        }

        auto bam_filenames = args[1 .. $];
        auto bam = new MultiBamReader(bam_filenames);
        enforce(bam.header.sorting_order == SortingOrder.coordinate,
                "All files must be coordinate-sorted");
        enforce(bam.has_index, "All files must be indexed");

        printer.bam = bam;

        string[string] rg2sm;
        foreach (rg; bam.header.read_groups) {
            rg2sm[rg.identifier] = rg.sample;
        }

        InputRange!(CustomBamRead) reads;
        if (bed_filename !is null) {
            auto bed = parseBed(bed_filename, bam);
	    bool simplify = mode == Mode.base;
            printer.setBed(parseBed(bed_filename, bam, simplify,
				    &printer.raw_bed_lines));
            reads = inputRangeObject(bam.getReadsOverlapping(bed)
                                     .map!(r => CustomBamRead(r, rg2sm)));
        } else {
            reads = inputRangeObject(bam.reads.map!(r => CustomBamRead(r, rg2sm)));
        }

        auto filtered_reads = inputRangeObject(filtered(reads, read_filter));
        auto pileup = pileupColumns(filtered_reads);

        int last_ref_id = -2;

        foreach (column; pileup) {
            auto ref_name = bam.reference_sequences[column.ref_id].name;

            if (column.ref_id != last_ref_id) {
                last_ref_id = column.ref_id;
                stderr.writeln("Processing reference #", column.ref_id + 1,
                               " (", ref_name, ")");
            }

            printer.push(column);
        }

        printer.close();
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
