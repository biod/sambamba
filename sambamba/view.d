/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
import bio.sam.reader;
import bio.core.region;

import sambamba.utils.common.filtering;
import sambamba.utils.common.progressbar;
import sambamba.utils.view.alignmentrangeprocessor;
import sambamba.utils.view.headerserializer;

import bio.core.utils.format;

import std.stdio;
import std.string;
import std.array;
import std.traits;
import std.getopt;
import std.parallelism;
import std.algorithm;
import std.math : isNaN;
import std.random;
import std.range;
import std.conv;

void printUsage() {
    stderr.writeln("Usage: sambamba-view [options] <input.bam | input.sam> [region1 [...]]");
    stderr.writeln();
    stderr.writeln("Options: -F, --filter=FILTER");
    stderr.writeln("                    set custom filter for alignments");
    stderr.writeln("         -f, --format=sam|bam|json");
    stderr.writeln("                    specify which format to use for output (default is SAM)");
    stderr.writeln("         -h, --with-header");
    stderr.writeln("                    print header before reads (always done for BAM output)");
    stderr.writeln("         -H, --header");
    stderr.writeln("                    output only header to stdout (if format=bam, the header is printed as SAM)");
    stderr.writeln("         -I, --reference-info");
    stderr.writeln("                    output to stdout only reference names and lengths in JSON");
    stderr.writeln("         -c, --count");
    stderr.writeln("                    output to stdout only count of matching records, hHI are ignored");
    stderr.writeln("         -v, --valid");
    stderr.writeln("                    output only valid alignments");
    stderr.writeln("         -S, --sam-input");
    stderr.writeln("                    specify that input is in SAM format");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("                    show progressbar in STDERR (works only for BAM files with no regions specified)");
    stderr.writeln("         -l, --compression-level");
    stderr.writeln("                    specify compression level (from 0 to 9, works only for BAM output)");
    stderr.writeln("         -o, --output-filename");
    stderr.writeln("                    specify output filename");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -s, --subsample=FRACTION");
    stderr.writeln("                    subsample reads (read pairs)");
    stderr.writeln("         --subsampling-seed=SEED");
    stderr.writeln("                    set seed for subsampling");
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
bool with_header;
bool header_only;
bool reference_info_only;
bool count_only;
bool skip_invalid_alignments;
bool is_sam;

bool show_progress;

int compression_level = -1;
string output_filename;
uint n_threads;
double subsample_frac; // NaN by default
uint subsample_threshold;

uint subsampling_seed;

// implementation of subsampling is the same as in samtools
uint simpleHash(string s) {
    uint h = 0;
    foreach (char c; s)
        h = (h << 5) - h + c;
    return h + subsampling_seed;
}

bool keepRead(R)(auto ref R read) {
    auto h = simpleHash(read.name);
    return (h&0xFFFFFF) < subsample_threshold;
}

InputRange!(ElementType!R) filteredView(R)(R reads) {
    if (isNaN(subsample_frac))
        return filtered(reads).inputRangeObject();

    return filtered(reads).filter!keepRead()
                          .inputRangeObject();
}

version(standalone) {
    int main(string[] args) {
        return view_main(args);
    }
}

int view_main(string[] args) {
    n_threads = totalCPUs - 1;

    subsampling_seed = unpredictableSeed;

    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "filter|F",            &query,
               "format|f",            &format,
               "with-header|h",       &with_header,
               "header|H",            &header_only,
               "reference-info|I",    &reference_info_only,
               "count|c",             &count_only,
               "valid|v",             &skip_invalid_alignments,
               "sam-input|S",         &is_sam,
               "show-progress|p",     &show_progress,
               "compression-level|l", &compression_level,
               "output-filename|o",   &output_filename,
               "nthreads|t",          &n_threads,
               "subsample|s",         &subsample_frac,
               "subsampling-seed",    &subsampling_seed);
        
        if (args.length < 2) {
            printUsage();
            return 0;
        }

        if (!isNaN(subsample_frac))
            subsample_threshold = (0x1000000 * subsample_frac).to!uint;

        auto task_pool = new TaskPool(n_threads);
        scope(exit) task_pool.finish();
        if (!is_sam) {
            auto bam = new BamReader(args[1], task_pool); 
            return sambambaMain(bam, task_pool, args);
        } else {
            auto sam = new SamReader(args[1]);
            return sambambaMain(sam, task_pool, args);
        }
    } catch (Exception e) {
        stderr.writeln("sambamba-view: ", e.msg);

        version(development) {
            throw e; // rethrow to see detailed message
        }

        return 1;
    }
}

mixin GlobalFilter;

File output_file() @property {
    if (_f == File.init) {
        if (output_filename is null)
            _f = stdout;
        else
            _f = File(output_filename, "w+");
    }
    return _f;
}

File _f;

// In fact, $(D bam) is either BAM or SAM file
int sambambaMain(T)(T _bam, TaskPool pool, string[] args) 
    if (is(T == SamReader) || is(T == BamReader)) 
{
    auto bam = _bam; // FIXME: uhm, that was a workaround for some closure-related bug

    if (reference_info_only && !count_only) {
        outputReferenceInfoJson(bam);
        return 0;
    }

    if (header_only && !count_only) {
        // write header to stdout
        (new HeaderSerializer(stdout, format)).writeln(bam.header);
    } else if (with_header && !count_only && format != "bam") {
        // for BAM, header will be written by writeBAM function
        (new HeaderSerializer(output_file, format)).writeln(bam.header);
    }
    
    if (header_only) {
        output_file.close();
        return 0;
    }

    read_filter = new NullFilter();

    if (skip_invalid_alignments) {
        read_filter = new AndFilter(read_filter, new ValidAlignmentFilter());
    }

    if (query !is null) {
        auto query_filter = createFilterFromQuery(query); 
        if (query_filter is null)
            return 1;
        read_filter = new AndFilter(read_filter, query_filter);
    }

    int processAlignments(AlignmentRangeProcessor)(AlignmentRangeProcessor processor) {
        static if (is(T == SamReader)) {
            if (args.length > 2) {
                stderr.writeln("sorry, accessing regions is unavailable for SAM input");
                return 1;
            }
        }

        if (args.length == 2) {
            
            static if (is(T == BamReader)) {
                if (show_progress) {
                    auto bar = new shared(ProgressBar)();
                    auto reads = bam.readsWithProgress((lazy float p) { bar.update(p); });
                    processor.process(filteredView(reads), bam);
                    bar.finish();
                } else {
                    if (cast(NullFilter) read_filter && isNaN(subsample_frac))
                        processor.process(bam.reads!withoutOffsets(), bam);
                    else
                        processor.process(filteredView(bam.reads!withoutOffsets), bam);
                }
            } else { // SamFile
                if (cast(NullFilter) read_filter && isNaN(subsample_frac))
                    processor.process(bam.reads, bam);
                else
                    processor.process(filteredView(bam.reads), bam);
            }
        } 

        // for BAM, random access is available
        static if (is(T == BamReader)) {
            if (args.length > 2) {
                auto regions = map!parseRegion(args[2 .. $]);

                alias ReturnType!(ReferenceSequence.opSlice) AlignmentRange;
                auto alignment_ranges = new AlignmentRange[regions.length];

                size_t i = 0;
                foreach (ref r; regions) {
                    alignment_ranges[i++] = bam[r.reference][r.beg .. r.end];
                }

                auto reads = joiner(alignment_ranges);
                processor.process(filteredView(reads), bam);
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
        if (format == "bam")             // will close the file
            return processAlignments(new BamSerializer(output_file, compression_level, pool));

        scope (exit) output_file.close();
        switch (format) {
            case "sam":
                return processAlignments(new SamSerializer(output_file, pool));
            case "json":
                return processAlignments(new JsonSerializer(output_file, pool));
            case "msgpack":
                return processAlignments(new MsgpackSerializer(output_file));
            default:
                stderr.writeln("output format must be one of sam, bam, json");
                return 1;
        }
    }

    return 0;
}
