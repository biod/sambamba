/*
    This file is part of Sambamba.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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
import bio.core.region;
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
    stderr.writeln("Usage: sambamba-depth region|window|base [options] input.bam  [input2.bam [...]]");
    stderr.writeln();
    stderr.writeln("          All BAM files must be coordinate-sorted and indexed.");
    stderr.writeln();
    stderr.writeln("          The tool has three modes: base, region, and window,");
    stderr.writeln("          each name means per which unit to print the statistics.");
    stderr.writeln();
    stderr.writeln("Common options:");
    stderr.writeln("         -F, --filter=FILTER");
    stderr.writeln("                    set custom filter for alignments; the default value is");
    stderr.writeln("                    'mapping_quality > 0 and not duplicate and not failed_quality_control'");
    stderr.writeln("         -o, --output-file=FILENAME");
    stderr.writeln("                    output filename (by default /dev/stdout)");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -c, --min-coverage=MINCOVERAGE");
    stderr.writeln("                    minimum mean coverage for output (default: 0 for region/window, 1 for base)");
    stderr.writeln("         -C, --max-coverage=MAXCOVERAGE");
    stderr.writeln("                    maximum mean coverage for output");
    stderr.writeln("         -q, --min-base-quality=QUAL");
    stderr.writeln("                    don't count bases with lower base quality");
    stderr.writeln("         --combined");
    stderr.writeln("                    output combined statistics for all samples");
    stderr.writeln("         -a, --annotate");
    stderr.writeln("                    add additional column of y/n instead of");
    stderr.writeln("                    skipping records not satisfying the criteria");
    stderr.writeln("         -m, --fix-mate-overlaps");
    stderr.writeln("                    detect overlaps of mate reads and handle them on per-base basis");
    stderr.writeln("base subcommand options:");
    stderr.writeln("         -L, --regions=FILENAME|REGION");
    stderr.writeln("                    list or regions of interest or a single region in form chr:beg-end (optional)");
    stderr.writeln("         -z, --report-zero-coverage (DEPRECATED, use --min-coverage=0 instead)");
    stderr.writeln("                    don't skip zero coverage bases");
    stderr.writeln("region subcommand options:");
    stderr.writeln("         -L, --regions=FILENAME|REGION");
    stderr.writeln("                    list or regions of interest or a single region in form chr:beg-end (required)");
    stderr.writeln("         -T, --cov-threshold=COVTHRESHOLD");
    stderr.writeln("                    multiple thresholds can be provided,");
    stderr.writeln("                    for each one an extra column will be added,");
    stderr.writeln("                    the percentage of bases in the region");
    stderr.writeln("                    where coverage is more than this value");
    stderr.writeln("window subcommand options:");
    stderr.writeln("         -w, --window-size=WINDOWSIZE");
    stderr.writeln("                    breadth of the window, in bp (required)");
    stderr.writeln("         --overlap=OVERLAP");
    stderr.writeln("                    overlap of successive windows, in bp (default is 0)");
    stderr.writeln("         -T, --cov-threshold=COVTHRESHOLD");
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

        size_t start_index = 0;
        size_t end_index = start_index;

        intervalTreeNode[][int] intervals;
        foreach (i, line; bed) {
            auto node = new intervalTreeNode(line.start, line.end, i);
            intervals[line.ref_id] ~= node;
        }

        trees_.length = reduce!max(intervals.keys) + 1;
        foreach (ref_id, ivs; intervals)
            trees_[ref_id] = new intervalTree(ivs);
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

enum MateOverlapStatus : ubyte {
    none = 0,
    detected = 1,
    fixed = 2,
    past = 3
}

struct CustomBamRead {
    MultiBamRead!BamRead read;
    alias read this;

    this(MultiBamRead!BamRead read, uint[string] rg2id) {
        this.read = read;
        auto rg_value = read["RG"];
        if (rg2id.length > 0 && !rg_value.is_nothing) {
            auto rg_str = *(cast(string*)(&rg_value));
            auto p = rg_str in rg2id;
            if (!p)
              throw new Exception("error in read " ~ read.name ~
                                  ": read group " ~ rg_str ~ " is not present in the header");
            sample_id = *p;
        }

        // FNV hashing
        ulong h = 14695981039346656037UL;
        foreach (ubyte b; read.name) {
            h ^= b;
            h *= 1099511628211UL;
        }
        read_name_hash = h;
    }

    ulong read_name_hash;
    uint sample_id;
    MateOverlapStatus mate_overlap = MateOverlapStatus.none;

    CustomBamRead dup() @property const {
        CustomBamRead r = void;
        r.read = read.dup;
        r.sample_id = sample_id;
        r.read_name_hash = read_name_hash;
        r.mate_overlap = mate_overlap;
        return r;
    }
}

alias Column = PileupRange!(InputRange!CustomBamRead).Column;

abstract class ColumnPrinter {
    double min_cov = 0.0;
    double max_cov = 1e50;
    ubyte min_base_quality = 0;
    MultiBamReader bam;

    bool combined = false;
    bool annotate = false;
    bool fix_mate_overlaps = false;
    File output_file;

    string[] sample_names;

    BamRegion[] raw_bed;

    void setBed(BamRegion[] bed) {
        raw_bed = bed;
    }

    string[] raw_bed_lines;

    abstract void init(ref string[] args);
    abstract void push(ref Column);
    abstract void close();

    uint getSampleId(R)(auto ref R read) {
        if (combined || sample_names.length == 1)
            return 0;
        return read.sample_id;
    }

    string getSampleName(uint sample_id) {
        return sample_names[sample_id];
    }

    private {
        Tuple!(ulong, size_t)[] read_name_hashes_buf;
        Tuple!(size_t, size_t)[] overlapping_mate_positions_buf;
    }

    Tuple!(size_t, size_t)[] overlapping_mate_positions;

    void detectOverlappingMates(ref Column column) {
        if (!fix_mate_overlaps)
            return;

        size_t n = column.coverage;

        if (n == 0) {
            overlapping_mate_positions = [];
            return;
        }

        if (n > read_name_hashes_buf.length) {
            read_name_hashes_buf.length = n;
            overlapping_mate_positions_buf.length = n;
        }

        for (size_t i = 0; i < n; i++)
            read_name_hashes_buf[i] = tuple(column.reads[i].read_name_hash, i);

        sort!((a, b) => a[0] < b[0])(read_name_hashes_buf[0 .. n]);

        size_t n_overlaps = 0;

        for (size_t i = 0; i < n - 1; i++) {
            if (read_name_hashes_buf[i][0] != read_name_hashes_buf[i + 1][0]) {
                auto idx = read_name_hashes_buf[i][1];
                if (column.reads[idx].mate_overlap != MateOverlapStatus.none)
                    column.reads[idx].mate_overlap = MateOverlapStatus.past;
                continue;
            }

            auto idx1 = read_name_hashes_buf[i][1];
            auto idx2 = read_name_hashes_buf[i + 1][1];

            if (column.reads[idx1].sample_id == column.reads[idx2].sample_id &&
                column.reads[idx1].name == column.reads[idx2].name)
            {
                if (column.reads[idx1].mate_overlap != MateOverlapStatus.none &&
                    column.reads[idx2].mate_overlap != MateOverlapStatus.none)
                {
                    assert(column.reads[idx1].mate_overlap == column.reads[idx2].mate_overlap);
                }

                overlapping_mate_positions_buf[n_overlaps++] = tuple(idx1, idx2);

                // don't touch status if it's already set
                if (column.reads[idx1].mate_overlap == MateOverlapStatus.none)
                    column.reads[idx1].mate_overlap = MateOverlapStatus.detected;

                if (column.reads[idx2].mate_overlap == MateOverlapStatus.none)
                    column.reads[idx2].mate_overlap = MateOverlapStatus.detected;

                // don't consider rare cases of >= 3 reads with the same name
                i += 1;
            }
        }

        auto idx = read_name_hashes_buf[n - 1][1];
        if (column.reads[idx].mate_overlap != MateOverlapStatus.none &&
            ((n == 1) ||
             (read_name_hashes_buf[n - 2][0] != read_name_hashes_buf[n - 1][0])))
            column.reads[idx].mate_overlap = MateOverlapStatus.past;

        // store positions of overlapping mates for this column
        overlapping_mate_positions = overlapping_mate_positions_buf[0 .. n_overlaps];
    }

    // selects a mate of better quality out of two overlapping candidates
    auto ref R selectBetterMate(R)(auto ref R m1, auto ref R m2) {
        if (m1.current_base == '-' || m2.current_base == '-') {
            // if either one hits an indel, look at the mapping quality
            return m1.mapping_quality > m2.mapping_quality ?  m1 : m2;
        } else {
            // otherwise, choose one with higher base quality
            return m1.current_base_quality > m2.current_base_quality ? m1 : m2;
        }
    }
}

final class PerBasePrinter : ColumnPrinter {
    NonOverlappingRegionStatsCollector stats_collector;
    bool report_zero_coverage;

    private {
        int _prev_ref_id = -2;
        long _prev_position;
        bool _bed_is_provided;
    }

    override void init(ref string[] args) {
        getopt(args,
                std.getopt.config.caseSensitive,
                "report-zero-coverage|z", &report_zero_coverage);
        if (report_zero_coverage)
            min_cov = 0;
        if (min_cov == 0)
            report_zero_coverage = true;

        output_file.write("REF\tPOS\tCOV\tA\tC\tG\tT\tDEL\tREFSKIP");
        if (!combined)
            output_file.write("\tSAMPLE");
        if (annotate)
            output_file.write("\tFLAG");
        output_file.writeln();
    }

    override void setBed(BamRegion[] bed) {
        super.setBed(bed);
        _bed_is_provided = true;
        stats_collector = new NonOverlappingRegionStatsCollector(bed);
    }

    private string[] tails;
    private void initTails() {
        if (!tails.empty) return;

        if (combined) {
            tails ~= "\t0\t0\t0\t0\t0\t0\t0";
            if (annotate)
                tails[0] = tails[0] ~ (min_cov > 0 ? "\tn" : "\ty");
        } else {
            foreach (sample_name; sample_names) {
                tails ~= "\t0\t0\t0\t0\t0\t0\t0\t" ~ sample_name;
                if (annotate)
                    tails.back = tails.back ~ (min_cov > 0 ? "\tn" : "\ty");
            }
        }
    }

    private void writeEmptyColumns(long ref_id, long start, long end) {
        if (min_cov > 0 && !annotate)
            return;
        auto ref_name = bam.reference_sequences[cast(uint)ref_id].name;

        initTails();

        if (!_bed_is_provided) {
            foreach (pos; start .. end) {
                foreach (tail; tails)
                    output_file.writeln(ref_name, '\t', pos, tail);
            }
        } else {
            if (raw_bed.empty || raw_bed.front.ref_id > cast(uint)ref_id)
                return;
            while (!raw_bed.empty && raw_bed.front.ref_id < cast(uint)ref_id)
                raw_bed.popFront();
            while (!raw_bed.empty && raw_bed.front.ref_id == ref_id) {
                if (raw_bed.front.fullyLeftOf(cast(uint)ref_id, cast(uint)start)) {
                    raw_bed.popFront();
                    continue;
                }
                auto from = max(start, raw_bed.front.start);
                auto to = min(end, raw_bed.front.end);
                if (from >= to)
                    break;
                foreach (pos; from .. to)
                    foreach (tail; tails)
                    output_file.writeln(ref_name, '\t', pos, tail);
                raw_bed.front.start = cast(uint)to;
                if (raw_bed.front.start >= raw_bed.front.end)
                    raw_bed.popFront();
            }
            stats_collector = new NonOverlappingRegionStatsCollector(raw_bed);
        }
    }

    private {
        size_t[] coverage;
        size_t[] deletions;
        size_t[] ref_skips;
    }

    private void writeColumn(ref Column c) {
        if (coverage.empty) {
            deletions.length = max(1, combined ? 1 : sample_names.length);
            ref_skips.length = deletions.length;
            coverage.length = 5 * deletions.length;
        }

        coverage[] = 0;
        deletions[] = 0;
        ref_skips[] = 0;

        void processCurrentBase(R)(auto ref R read) {
            auto sample_id = getSampleId(read);
            if (read.current_base == '-') {
                if (read.cigar_operation.type == 'D')
                    deletions[sample_id] += 1;
                else
                    ref_skips[sample_id] += 1;
                return;
            }

            if (read.current_base_quality >= min_base_quality)
                coverage[5 * sample_id + Base5(read.current_base).internal_code] += 1;
        }

        detectOverlappingMates(c);

        foreach (read; c.reads) {
            if (read.mate_overlap == MateOverlapStatus.detected)
                continue; // process overlapping mates separately
            else
                processCurrentBase(read);
        }

        foreach (pair; overlapping_mate_positions) {
            auto read = selectBetterMate(c.reads[pair[0]], c.reads[pair[1]]);
            processCurrentBase(read);
        }

        foreach (sample_id; 0 .. coverage.length / 5) {
            auto cov = coverage[sample_id * 5 .. $][0 .. 5];
            size_t total_coverage = reduce!`a+b`(cov[]) +
                deletions[sample_id] + ref_skips[sample_id];

            bool ok = total_coverage >= min_cov && total_coverage <= max_cov;
            if (!ok && !annotate)
                return;

            output_file.write(bam.reference_sequences[c.ref_id].name, '\t',
                    c.position, '\t', total_coverage);
            foreach (i; 0 .. 4)
                output_file.write('\t', cov[i]);
            output_file.write('\t', deletions[sample_id], '\t', ref_skips[sample_id]);

            if (!combined)
                output_file.write('\t', getSampleName(sample_id.to!uint));

            if (annotate)
                output_file.write('\t', ok ? 'y' : 'n');
            output_file.writeln();
        }
    }

    private bool outputRequired(int ref_id, ulong position) {
        if (stats_collector is null)
            return true;
        bool output = false;
        stats_collector.nextColumn(cast(uint)ref_id, cast(uint)position,
                                   (size_t id) { output = true; });
        return output;
    }

    override void push(ref Column c) {
        if (min_cov > 0) {
            if (outputRequired(c.ref_id, c.position))
                writeColumn(c);
            return;
        }

        if (_prev_ref_id == -2) {
            foreach (id; 0 .. c.ref_id)
                writeEmptyColumns(id, 0, bam.reference_sequences[id].length);
            writeEmptyColumns(c.ref_id, 0, c.position);
        } else if (_prev_ref_id != c.ref_id) {
            writeEmptyColumns(_prev_ref_id, _prev_position + 1,
                    bam.reference_sequences[_prev_ref_id].length);
            writeEmptyColumns(c.ref_id, 0, c.position);
        } else if (_prev_position != c.position - 1) {
            writeEmptyColumns(c.ref_id, _prev_position + 1, c.position);
        }

        _prev_ref_id = c.ref_id;
        _prev_position = c.position;

        if (outputRequired(c.ref_id, c.position))
            writeColumn(c);
    }

    override void close() {
        if (!report_zero_coverage)
            return;

        if (_prev_ref_id == -2) {
            foreach (id; 0 .. bam.reference_sequences.length)
                writeEmptyColumns(id, 0, bam.reference_sequences[id].length);
        } else {
            writeEmptyColumns(_prev_ref_id, _prev_position + 1,
                              bam.reference_sequences[_prev_ref_id].length);
            foreach (id; _prev_ref_id + 1 .. bam.reference_sequences.length)
                writeEmptyColumns(id, 0, bam.reference_sequences[id].length);
        }
    }
}

final class PerSampleRegionData {
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

abstract class PerRegionPrinter : ColumnPrinter {
    RegionStatsCollector stats_collector;
    private PerSampleRegionData[] samples;

    static immutable default_bed_fields = ["chrom", "chromStart", "chromEnd"];

    private void printBedHeader(size_t n_before) {
        output_file.write("# ");
        foreach (field; default_bed_fields[0 .. min($, n_before)])
            output_file.write(field, "\t");
        foreach (k; 3 .. n_before)
            output_file.write("F", k, "\t");
        output_file.write("readCount\tmeanCoverage");
        foreach (cov; cov_thresholds)
            output_file.write("\tpercentage", cov);
        if (!combined)
            output_file.write("\tsampleName");

        if (annotate)
            output_file.write("\tmeanCovWithinBounds");
        output_file.write("\n");
        output_file.flush();
    }

    private void countRead(R)(auto ref R read, size_t id) {
        auto n = countOverlappingBases(read, id);
        auto data = getSampleData(getSampleId(read));
        data.n_bases(id) += n;

        // count the read only if at least one base is good
        if (n > 0)
            data.n_reads(id) += 1;
    }

    private size_t countOverlappingBases(R)(auto ref R read, size_t id, ulong start_pos=0) {
        auto sample_id = getSampleId(read);
        assert(sample_names.empty || sample_id < sample_names.length, "Invalid sample ID");
        auto data = getSampleData(sample_id);

        auto region = getRegionById(id);

        auto pos = read.position; // current position on the reference
        auto q = read.base_qualities;
        size_t n; // number of read bases that are not insertions
                  // and also have good quality
        foreach (op; read.cigar) {
            if (op.is_match_or_mismatch)
                foreach (qual; q[0 .. min(op.length, $)]) {
                    n += region.overlaps(region.ref_id, pos) &&
                         qual >= min_base_quality && pos >= start_pos;
                    ++pos;
                }
            else if (op.is_reference_consuming)
                pos += op.length;

            // min(op.length, $) protects from invalid memory accesses
            // possible only when the input data is incorrect
            if (op.is_query_consuming)
                q = q[min(op.length, $) .. $];
        }
        return n;
    }

    abstract BamRegion getRegionById(size_t id);
    abstract PerSampleRegionData getSampleData(uint sample_id);
    abstract bool isFirstOccurrence(size_t id);
    abstract void markAsSeen(size_t id);
    abstract void writeOriginalBedLine(size_t id);

    uint[] cov_thresholds;
    bool[] is_first_occurrence;

    uint[] cov_per_sample;

    override void init(ref string[] args) {
        getopt(args,
               std.getopt.config.caseSensitive,
               "cov-threshold|T", &cov_thresholds);
    }

    private void uncountOverlappingMates(R)(ref R r1, ref R r2, size_t id, ulong curr_pos) {
        if (r1.mate_overlap == MateOverlapStatus.fixed &&
            r2.mate_overlap == MateOverlapStatus.fixed)
            return;

        assert(r1.mate_overlap == MateOverlapStatus.detected);
        assert(r2.mate_overlap == MateOverlapStatus.detected);

        // re-count all good bases
        auto n1_full = countOverlappingBases(r1, id);
        auto n2_full = countOverlappingBases(r2, id);

        // now count bases starting from the current column
        auto n1 = r1.position == curr_pos ? n1_full : countOverlappingBases(r1, id, curr_pos);
        auto n2 = r2.position == curr_pos ? n2_full : countOverlappingBases(r2, id, curr_pos);

        auto data = getSampleData(r1.sample_id);

        data.n_bases(id) -= n1 + n2; // this count is then dealt with on per-base basis

        data.n_reads(id) -= (n1_full > 0) + (n2_full > 0);
        data.n_reads(id) += (n1_full + n2_full > 0); // count only one read instead of two

        // don't set status to fixed just yet, there might be other regions as well
    }

    private void fixRegionBaseCounter(ref Column column, size_t region_id) {
        foreach (p; overlapping_mate_positions)
            uncountOverlappingMates(column.reads[p[0]],
                                    column.reads[p[1]], region_id, column.position);
    }

    private void markOverlappingMatesAsFixed(ref Column column) {
        foreach (p; overlapping_mate_positions) {
            assert(column.reads[p[0]].mate_overlap != MateOverlapStatus.past);
            assert(column.reads[p[1]].mate_overlap != MateOverlapStatus.past);
            column.reads[p[0]].mate_overlap = MateOverlapStatus.fixed;
            column.reads[p[1]].mate_overlap = MateOverlapStatus.fixed;
        }
    }

    override void push(ref Column column) {
        uint ref_id = column.ref_id.to!uint;
        pos_t position = column.position.to!pos_t;

        if (cov_per_sample.length == 0) {
            cov_per_sample.length = max(1, combined ? 1 : sample_names.length);
        }

        detectOverlappingMates(column);

        void processCurrentBase(R)(auto ref R read, size_t region_id) {
            if (read.current_base_quality < min_base_quality)
                return;
            auto sample_id = getSampleId(read);
            auto data = getSampleData(sample_id);
            data.n_bases(region_id) += 1;
            cov_per_sample[sample_id] += 1;
        }

        void countPreviouslySeenMateOverlaps(size_t id) {
            foreach (pair; overlapping_mate_positions) {
                auto m1 = column.reads[pair[0]];
                auto m2 = column.reads[pair[1]];
                if (m1.mate_overlap != MateOverlapStatus.fixed)
                    continue;

                // don't set n_bases here, only n_reads
                auto n1 = countOverlappingBases(m1, id);
                auto n2 = countOverlappingBases(m2, id);

                if (n1 + n2 == 0)
                    continue;

                auto data = getSampleData(m1.sample_id);
                data.n_reads(id) += 1;
            }
        }

        bool fixes_applied = false;

        stats_collector.nextColumn(ref_id, position,
            (size_t id) {
                if (isFirstOccurrence(id)) {
                    foreach (read; column.reads)
                        if (read.mate_overlap != MateOverlapStatus.fixed)
                            countRead(read, id);
                    countPreviouslySeenMateOverlaps(id);
                    markAsSeen(id);
                } else {
                    foreach (read; column.reads_starting_here)
                        countRead(read, id);
                }

                fixRegionBaseCounter(column, id);
                fixes_applied = true;

                cov_per_sample[] = 0;

                foreach (ref read; column.reads) {
                    if (read.mate_overlap != MateOverlapStatus.none)
                    {
                        if (read.mate_overlap != MateOverlapStatus.past)
                            continue;
                        processCurrentBase(read, id);
                    } else {
                        if (read.current_base_quality >= min_base_quality)
                            cov_per_sample[getSampleId(read)] += 1;
                    }
                }

                foreach (pair; overlapping_mate_positions) {
                    auto m1 = column.reads[pair[0]], m2 = column.reads[pair[1]];
                    processCurrentBase(selectBetterMate(m1, m2), id);
                }

                foreach (sample_id; iota(cov_per_sample.length.to!uint)) {
                    auto data = getSampleData(sample_id);
                    foreach (i, threshold; cov_thresholds)
                        if (cov_per_sample[sample_id] >= threshold)
                            data.coverage_count(i, id) += 1;
                }
        });

        if (fixes_applied)
            markOverlappingMatesAsFixed(column);
    }

    void printRegionStats(uint sample_id, size_t id, PerSampleRegionData data) {
        auto region = getRegionById(id);
        auto length = region.end - region.start;
        with(output_file) {
            auto mean_cov = data.n_bases(id).to!float / length;

            bool ok = mean_cov >= this.min_cov && mean_cov <= this.max_cov;

            if (!ok && !this.annotate)
                return;

            writeOriginalBedLine(id);
            write(data.n_reads(id), '\t', mean_cov);
            foreach (j; 0 .. cov_thresholds.length) {
                auto percentage = data.coverage_count(j, id).to!float * 100 / length;
                if (cov_thresholds[j] == 0)
                    percentage = 100.0;
                write('\t', percentage);
            }

            if (!combined)
                write('\t', getSampleName(sample_id));

            if (annotate)
                write('\t', ok ? 'y' : 'n');

            writeln();
            flush();
        }
    }
}

final class PerBedRegionPrinter : PerRegionPrinter {
    bool[] is_first_occurrence;

    override PerSampleRegionData getSampleData(uint id) {
        if (samples.length == 0) {
            samples.length = max(1, combined ? 1 : sample_names.length);
            foreach (k; 0 .. samples.length)
                samples[k] = new PerSampleRegionData(cov_thresholds.length, raw_bed.length);
        }
        assert(id < samples.length, "Invalid sample ID: " ~
               id.to!string ~ "/" ~ samples.length.to!string);
        assert(samples[id] !is null);
        return samples[id];
    }

    override bool isFirstOccurrence(size_t id) {
        return is_first_occurrence[id];
    }

    override void markAsSeen(size_t id) {
        is_first_occurrence[id] = false;
    }

    override void writeOriginalBedLine(size_t id) {
        import std.string;
        raw_bed_lines[id] = std.string.stripRight(raw_bed_lines[id]);
        output_file.write(raw_bed_lines[id], "\t");
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

        printBedHeader(raw_bed_lines[0].split().length);
    }

    override void close() {
        foreach (id; 0 .. raw_bed.length) {
            foreach (sample_id; iota(samples.length.to!uint))
                printRegionStats(sample_id, id, getSampleData(sample_id));
        }
    }
}

final class PerWindowPrinter : PerRegionPrinter {
    size_t window_size;
    size_t overlap;
    size_t step;
    size_t n; // number of windows to keep at each moment

    bool[] is_first_occurrence;

    size_t leftmost_window_index = 0;
    size_t leftmost_window_start_pos = 0;
    int window_ref_id = -1;
    size_t ref_length;

    void printWindowStats(size_t id) {
        foreach (sample_id; iota(samples.length.to!uint))
            printRegionStats(sample_id, leftmost_window_index, getSampleData(sample_id));
    }

    void resetAllWindows() {
        foreach (data; samples) {
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

    override PerSampleRegionData getSampleData(uint id) {
        if (samples.length == 0) {
            samples.length = max(1, combined ? 1 : sample_names.length);
            foreach (k; 0 .. samples.length)
                samples[k] = new PerSampleRegionData(cov_thresholds.length, n);
        }
        return samples[id];
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

    override void writeOriginalBedLine(size_t id) {
        auto region = getRegionById(id);
        output_file.write(bam.reference_sequences[region.ref_id].name, '\t',
                          region.start, '\t',
                          region.end, '\t');
    }

    override void init(ref string[] args) {
        getopt(args, std.getopt.config.caseSensitive,
               "window-size|w", &window_size,
               "overlap", &overlap,
               "cov-threshold|T", &cov_thresholds);

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

        printBedHeader(3);
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
            while (leftmost_window_start_pos + window_size <= ref_length)
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
        while (leftmost_window_start_pos + window_size <= ref_length)
            finishLeftMostWindow();
        foreach (k; window_ref_id + 1 .. bam.reference_sequences.length)
            printEmptyWindows(cast(int)k);
        output_file.close();
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

    if (mode == Mode.base)
        printer.min_cov = 1;

    string output_fn;

    args = args[1 .. $];

    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               std.getopt.config.passThrough,
               "filter|F",               &query,
               "output-filename|o",      &output_fn,
               "nthreads|t",             &n_threads,
               "min-coverage|c",         &printer.min_cov,
               "max-coverage|C",         &printer.max_cov,
               "min-base-quality|q",     &printer.min_base_quality,
               "annotate|a",             &printer.annotate,
               "combined",               &printer.combined,
               "fix-mate-overlaps|m",    &printer.fix_mate_overlaps);

        if (output_fn is null)
            printer.output_file = stdout;
        else
            printer.output_file = File(output_fn, "w+");

        if (mode != Mode.window) {
            getopt(args,
                   std.getopt.config.caseSensitive,
                   std.getopt.config.passThrough,
                   "regions|L", &bed_filename);
        }

        if (mode == Mode.region && bed_filename is null) {
            stderr.writeln("BED file or a region must be provided in region mode");
            return 1;
        }

        // handles subcommand arguments and removes them from the list
        printer.init(args);

        defaultPoolThreads = max(n_threads, 0);

        if (query !is null) {
            read_filter = createFilterFromQuery(query);
        } else {
            read_filter = createFilterFromQuery("mapping_quality > 0 and "
                    "not duplicate and "
                    "not failed_quality_control");
        }

        auto bam_filenames = args[1 .. $];
        auto bam = new MultiBamReader(bam_filenames);
        enforce(bam.header.sorting_order == SortingOrder.coordinate,
                "All files must be coordinate-sorted");
        enforce(bam.has_index, "All files must be indexed");

        printer.bam = bam;

        uint[string] sm2id;
        uint[string] rg2id;
        foreach (rg; bam.header.read_groups) {
            if (rg.sample !in sm2id) {
                sm2id[rg.sample] = printer.sample_names.length.to!uint;
                printer.sample_names ~= rg.sample;
            }
            rg2id[rg.identifier] = sm2id[rg.sample];
        }

        if (printer.sample_names.empty)
            printer.sample_names = ["*"];

        InputRange!(CustomBamRead) reads;
        if (bed_filename !is null) {
            BamRegion[] bed;
            try {
                bed = parseBed(bed_filename, bam);
                if (mode == Mode.base) {
                    parseBed(bed_filename, bam, true, &printer.raw_bed_lines);
                    printer.setBed(parseBed(bed_filename, bam, true));
                } else {
                    printer.setBed(parseBed(bed_filename, bam, false, &printer.raw_bed_lines));
                }
            } catch (Exception e) {
                auto region = parseRegion(bed_filename);
                enforce(bam.hasReference(region.reference),
                    "couldn't open file " ~ bed_filename ~
                    " or find reference " ~ region.reference);
                auto ref_id = bam[region.reference].id;
                bed ~= BamRegion(ref_id, region.beg, region.end);
                if (bed[0].end == uint.max)
                    bed[0].end = bam[region.reference].length;
                auto s = region.reference ~ "\t" ~
                         bed[0].start.to!string() ~ "\t" ~
                         bed[0].end.to!string();
                printer.raw_bed_lines = [s];
                printer.setBed(bed);
            }
            reads = inputRangeObject(bam.getReadsOverlapping(bed)
                    .map!(r => CustomBamRead(r, rg2id)));
        } else {
            reads = inputRangeObject(bam.reads.map!(r => CustomBamRead(r, rg2id)));
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
