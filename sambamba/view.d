/*
    This file is part of Sambamba.
    Copyright © 2012-2016    Artem Tarasov <lomereiter@gmail.com>
    Copyright © 2012-2017    Pjotr Prins <pjotr.prins@thebird.nl>

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
import bio.bam.read;
import bio.bam.region;
import bio.bam.splitter;

import bio.core.region;
import bio.core.utils.format;

import bio.sam.reader;
import cram.reader;

import sambamba.utils.common.filtering;
import sambamba.utils.common.overwrite;
import sambamba.utils.common.progressbar;
import sambamba.utils.view.alignmentrangeprocessor;
import sambamba.utils.view.headerserializer;
import sambamba.utils.common.bed;
import sambamba.bio2.unpack;
// import core.sys.posix.stdlib; // for exit

import utils.version_ : addPG;

import std.array;
import std.algorithm;
import std.conv;
import std.exception;
import std.getopt;
import std.parallelism;
import std.math : isNaN;
import std.random;
import std.range;
import std.stdio;
import std.string;
import std.traits;

void printUsage() {
    stderr.writeln("Usage: sambamba-view [options] <input.bam | input.sam> [region1 [...]]

Options: -F, --filter=FILTER
                    set custom filter for alignments
         --num-filter=NUMFILTER
                    filter flag bits; 'i1/i2' corresponds to -f i1 -F i2 samtools arguments;
                    either of the numbers can be omitted
         -f, --format=sam|bam|cram|json|unpack
                    specify which format to use for output (default is SAM);
                    unpack streams unpacked BAM
         -h, --with-header
                    print header before reads (always done for BAM output)
         -H, --header
                    output only header to stdout (if format=bam, the header is printed as SAM)
         -I, --reference-info
                    output to stdout only reference names and lengths in JSON
         -L, --regions=FILENAME
                    output only reads overlapping one of regions from the BED file
         -c, --count
                    output to stdout only count of matching records, hHI are ignored
         -v, --valid
                    output only valid alignments
         -S, --sam-input
                    specify that input is in SAM format
         -C, --cram-input
                    specify that input is in CRAM format
         -T, --ref-filename=FASTA
                    specify reference for writing CRAM
         -p, --show-progress
                    show progressbar in STDERR (works only for BAM files with no regions specified)
         -l, --compression-level
                    specify compression level (from 0 to 9, works only for BAM output)
         -o, --output-filename
                    specify output filename
         -t, --nthreads=NTHREADS
                    maximum number of threads to use
         -s, --subsample=FRACTION
                    subsample reads (read pairs)
         --subsampling-seed=SEED
                    set seed for subsampling");
}

void outputReferenceInfoJson(T)(T bam) {

    auto w = stdout.lockingTextWriter;
    w.put('[');

    bool first = true;
    foreach (refseq; bam.reference_sequences) {
        if (first) {
            first = false;
        } else {
            w.put(',');
        }
        w.put(`{"name":`);
        writeJson((const(char)[] s) { w.put(s); }, refseq.name);
        w.put(`,"length":`);
        writeJson((const(char)[] s) { w.put(s); }, refseq.length);
        w.put('}');
    }

    w.put("]\n");
}

string format = "sam";
string query;
string numfilter;
string ref_fn;
bool with_header;
bool header_only;
bool reference_info_only;
bool count_only;
bool skip_invalid_alignments;
bool is_sam;
bool is_cram;

bool show_progress;

int compression_level = -1;
string output_filename;
uint n_threads;
double subsample_frac; // NaN by default
ulong subsampling_seed;
string bed_filename;

string[] unparsed_args;

version(standalone) {
    int main(string[] args) {
        return view_main(args);
    }
}

int view_main(string[] args) {
    foreach(arg; args) {
      if (arg == "--throw-error") {
        // undocumented: can throw a null pointer exception for testing debugger(s)
        char *p = null;
        *p = 'X'; // force an exception
      }
    }

    n_threads = totalCPUs;

    subsampling_seed = unpredictableSeed;
    subsampling_seed <<= 32;
    subsampling_seed += unpredictableSeed;

    unparsed_args = args.dup;

    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "filter|F",            &query,
               "num-filter",          &numfilter,
               "format|f",            &format,
               "with-header|h",       &with_header,
               "header|H",            &header_only,
               "reference-info|I",    &reference_info_only,
               "regions|L",           &bed_filename,
               "count|c",             &count_only,
               "valid|v",             &skip_invalid_alignments,
               "sam-input|S",         &is_sam,
               "cram-input|C",        &is_cram, // TODO: autodetection of format
               "show-progress|p",     &show_progress,
               "compression-level|l", &compression_level,
               "output-filename|o",   &output_filename,
               "nthreads|t",          &n_threads,
               "subsample|s",         &subsample_frac,
               "subsampling-seed",    &subsampling_seed,
               "ref-filename|T",      &ref_fn);

        if (args.length < 2) {
            printUsage();
            return 0;
        }

        // <cludge> cram writer doesn't share the task pool
        if (!is_sam && format == "cram")
          n_threads /= 2;
        // </cludge>

        protectFromOverwrite(args[1], output_filename);

        File output_file;
        if (output_filename is null || output_filename == "-")
          output_file = stdout;
        else
          output_file = File(output_filename, "w+");

        if (is_cram && is_sam)
            throw new Exception("only one of --sam-input and --cram-input can be specified");

        if (format == "unpack" ) {
          enforce(!is_cram && !is_sam,"we can only unpack bam files at this point!");
          // unpack is part of sambamba2
          return unpack_bams(args[1..$],output_file); // exits
        }
        auto task_pool = new TaskPool(n_threads);
        scope(exit) task_pool.finish();
        if (is_sam) {
          auto sam = new SamReader(args[1]);
          return sambambaMain(sam, task_pool, args, output_file);
        } else if (!is_cram) {
          auto bam = new BamReader(args[1], task_pool);
          return sambambaMain(bam, task_pool, args, output_file);
        } else {
          auto cram = new CramReader(args[1], task_pool);
          return sambambaMain(cram, task_pool, args, output_file);
        }
    } catch (Exception e) {
          stderr.writeln("sambamba-view: ", e.msg);

          version(development) {
            throw e; // rethrow to see detailed message
          }

          return 1;
        }
    }

auto filtered(R)(R reads, Filter f) {
  return reads.zip(f.repeat()).filter!q{a[1].accepts(a[0])}.map!q{a[0]}();
}

 int sambambaMain(T)(T _bam, TaskPool pool, string[] args, File output_file)
 {
    auto bam = _bam; // FIXME: uhm, that was a workaround for some closure-related bug

    if (reference_info_only && !count_only) {
        outputReferenceInfoJson(bam);
        return 0;
    }

    if (!header_only)  // = some processing is done
        addPG("view", unparsed_args, bam.header);

    if (header_only && !count_only) {
        // write header to stdout
        (new HeaderSerializer(stdout, format)).writeln(bam.header);
    } else if (with_header && !count_only &&
               format != "bam" && format != "cram")
    {
        (new HeaderSerializer(output_file, format)).writeln(bam.header);
    }

    if (header_only) {
        output_file.close();
        return 0;
    }

    Filter read_filter = new NullFilter();

    if (skip_invalid_alignments) {
        read_filter = new AndFilter(read_filter, new ValidAlignmentFilter());
    }

    if (numfilter !is null) {
        ushort i1, i2;
        auto masks = numfilter.splitter("/").array();
        if (masks.length > 0 && masks[0].length > 0) i1 = masks[0].to!ushort;
        if (masks.length > 1 && masks[1].length > 0) i2 = masks[1].to!ushort;
        read_filter = new AndFilter(read_filter, new FlagBitFilter(i1, i2));
    }

    if (query !is null) {
        auto query_filter = createFilterFromQuery(query);
        if (query_filter is null)
            return 1;
        read_filter = new AndFilter(read_filter, query_filter);
    }

    if (!isNaN(subsample_frac)) {
        auto subsample_filter = new SubsampleFilter(subsample_frac, subsampling_seed);
        read_filter = new AndFilter(subsample_filter, read_filter);
    }

    int processAlignments(P)(P processor) {
        static if (is(T == SamReader)) {
            if (args.length > 2) {
                stderr.writeln("region queries are unavailable for SAM input");
                return 1;
            }
        }

        void runProcessor(SB, R, F)(SB bam, R reads, F filter) {
            if (processor.is_serial)
                bam.assumeSequentialProcessing();
            if (cast(NullFilter) filter)
                processor.process(reads, bam);
            else
                processor.process(reads.filtered(filter), bam);
        }

        bool output_all_reads = args.length == 2 &&
          (bed_filename.empty || bam.header.sorting_order != SortingOrder.coordinate);
        static if (is(T == SamReader))
          output_all_reads = true;

        if (bed_filename.length > 0 && args.length > 2) {
            throw new Exception("specifying both region and BED filename is disallowed");
        }

        if (output_all_reads) {
            if (bed_filename !is null) {
                auto regions = parseBed(bed_filename, bam);
                read_filter = new AndFilter(read_filter, new BedFilter(regions));
            }

            static if (is(T == BamReader)) {
                if (show_progress) {
                    auto bar = new shared(ProgressBar)();
                    auto reads = bam.readsWithProgress((lazy float p) { bar.update(p); });
                    runProcessor(bam, reads, read_filter);
                    bar.finish();
                } else {
                    runProcessor(bam, bam.reads!withoutOffsets(), read_filter);
                }
            } else { // SamFile
                runProcessor(bam, bam.reads, read_filter);
            }
        } else {
        // for BAM/CRAM, random access is available
        static if (is(T == BamReader) || is(T == CramReader)) {
            if (args.length > 2) {
                auto regions = map!parseRegion(args[2 .. $]);

                alias typeof(bam.unmappedReads().front) Read;
                alias InputRange!Read AlignmentRange;
                auto alignment_ranges = new AlignmentRange[regions.length];

                size_t i = 0;
                foreach (region_description; args[2 .. $]) {
                    AlignmentRange range;
                    if (region_description == "*") {
                        range = bam.unmappedReads().inputRangeObject;
                    } else {
                        auto r = parseRegion(region_description);
                        auto seq = bam[r.reference];
                        if (r.end == uint.max)
                            r.end = seq.length;
                        range = seq[r.beg .. r.end].inputRangeObject;
                    }
                    alignment_ranges[i++] = range;
                }

                auto reads = joiner(alignment_ranges);
                runProcessor(bam, reads, read_filter);
            } else if (bed_filename.length > 0) {
                static if (is(T == CramReader)) {
                    throw new Exception("BED support is unavailable for CRAM");
                } else {
                    auto regions = parseBed(bed_filename, bam);
                    auto reads = bam.getReadsOverlapping(regions);
                    runProcessor(bam, reads, read_filter);
                }
            }
        }
        }

        return 0;
    }

    if (count_only) {
        auto counter = new ReadCounter();

        if (processAlignments(counter))
            return 1;
        writeln(counter.number_of_reads);
    } else {
        if (format == "bam") {
            // FIXME: dirty hack to avoid closing file twice, depends on std.stdio.File internals
            scope(exit) *(*cast(void***)(&output_file)) = null;
            return processAlignments(new BamSerializer(output_filename,
                                                       output_file,
                                                       compression_level, pool));
        } else if (format == "cram") {
            output_file.close();
            return processAlignments(new CramSerializer(output_filename, ref_fn,
                                                        n_threads));
        }

        scope (exit) output_file.close();
        switch (format) {
            case "sam":
                return processAlignments(new SamSerializer(output_file, pool));
            case "json":
                return processAlignments(new JsonSerializer(output_file, pool));
            case "msgpack":
                return processAlignments(new MsgpackSerializer(output_file));
            default:
                stderr.writeln("output format must be one of sam, bam, cram, json");
                return 1;
        }
    }

    return 0;
}
