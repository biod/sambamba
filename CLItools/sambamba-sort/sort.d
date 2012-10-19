/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.sort;

import bamfile;
import bamoutput;
import samheader;
import alignment;
import splitter;
import utils.tmpfile;

import std.range;
import std.algorithm;
import std.traits;
import std.array;
import std.ascii;
import std.numeric;
import std.parallelism;
import std.getopt;
import std.path;
import std.file;
import std.stream;
import std.stdio;
import core.atomic;

import common.comparators;
import common.nwayunion : nWayUnion;
import common.progressbar;

import thirdparty.mergesort;

void printUsage() {
    stderr.writeln("Usage: sambamba-sort [options] <input.bam>");
    stderr.writeln();
    stderr.writeln("Options: -m, --memory-limit=LIMIT");
    stderr.writeln("               approximate memory limit (it's not guaranteed that it won't be exceeded, because of garbage collection)");
    stderr.writeln("         --tmpdir=TMPDIR");
    stderr.writeln("               directory for storing intermediate files; default is system directory for temporary files");
    stderr.writeln("         -o, --out=OUTPUTFILE");
    stderr.writeln("               output file name; if not provided, the result is written to a file with .sorted.bam extension");
    stderr.writeln("         -n, --sort-by-name");
    stderr.writeln("               sort by read name instead of coordinate");
    stderr.writeln("         -l, --compression-level=COMPRESSION_LEVEL");
    stderr.writeln("               level of compression for sorted BAM, from 0 to 9");
    stderr.writeln("         -u, --uncompressed-chunks");
    stderr.writeln("               write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("               show progressbar in STDERR");
}

version(standalone) {
    int main(string[] args) {
        return sort_main(args);
    }
}

private bool sort_by_name;
private bool show_progress;

private shared(ProgressBar) bar;
private shared(float[]) weights;
private shared(float[]) merging_progress;

int sort_main(string[] args) {

    if (args.length < 2) {
        printUsage();
        return 0;
    }

    try {
        string memory_limit_str = null;

        size_t memory_limit = 512 * 1024 * 1024;
        string tmpdir = null;
        bool uncompressed_chunks;
        int compression_level = -1;

        string output_filename;

        getopt(args,
               std.getopt.config.caseSensitive,
               "memory-limit|m",        &memory_limit_str,
               "tmpdir",                &tmpdir,
               "out|o",                 &output_filename,
               "sort-by-name|n",        &sort_by_name,
               "uncompressed-chunks|u", &uncompressed_chunks,
               "compression-level|l",   &compression_level,
               "show-progress|p",       &show_progress);

        if (output_filename is null) {
            output_filename = setExtension(args[1], "sorted.bam");
        }

        if (memory_limit_str !is null) {
            memory_limit = parseMemory(memory_limit_str);
        }

        // TODO: find a better formula
        memory_limit /= 4;

        auto task_pool = new TaskPool(totalCPUs);
        scope(exit) task_pool.finish();
        auto bam = BamFile(args[1], task_pool);
       
        // ------------------ set sorting order in the header ----------------------

        auto header = bam.header;
        header.sorting_order = sort_by_name ? SortingOrder.queryname :
                                              SortingOrder.coordinate;
        auto header_text = toSam(header);

        // ----------------------- sort chunks -------------------------------------

        static Alignment[] sortChunk(Alignment[] chunk) {
            if (!sort_by_name) {
                mergeSort!compareAlignmentCoordinates(chunk, true); // threaded
            } else {
                mergeSort!compareReadNames(chunk, true);
            }
            return chunk;
        }

        string[] tmpfiles;
        auto num_of_chunks = 0;

        void writeSortedChunks(R)(R alignments) {
            foreach (chunk; map!sortChunk(chunks(alignments, memory_limit)))
            {
                auto fn = tmpFile(chunkBaseName(args[1], num_of_chunks), tmpdir);
                tmpfiles ~= fn;

                Stream stream = new BufferedFile(fn, FileMode.Out);
                scope(exit) stream.close();

                writeBAM(stream, header_text, bam.reference_sequences, 
                         chunk, uncompressed_chunks ? 0 : 1, task_pool);

                num_of_chunks += 1;
            }
        }

        if (show_progress) {
            stderr.writeln("Writing sorted chunks to temporary directory...");
            bar = new shared(ProgressBar)();
            auto alignments = bam.alignmentsWithProgress(
                                  (lazy float p){ bar.update(p); }
                              );
            writeSortedChunks(alignments);
            bar.finish();
        } else {
            writeSortedChunks(bam.alignments!withoutOffsets);
        }

        scope(exit) {
            // ---------------- remove temporary files at exit ---------------------
            foreach (tmpfile; tmpfiles) {
                remove(tmpfile);
            }
        }

        // ---------------------- merge sorted chunks ------------------------------

        // half of memory is for input buffers
        // and another half is for output buffers
        Stream stream = new BufferedFile(output_filename, FileMode.OutNew, 
                                         memory_limit / 2);
        scope(exit) stream.close();

        if (show_progress) {
            stderr.writeln("Merging sorted chunks...");
            weights.length = num_of_chunks;
            merging_progress.length = num_of_chunks;
            merging_progress[] = 0.0;

            alias ReturnType!(BamFile.alignmentsWithProgress!withoutOffsets) AlignmentRangePB;
            auto alignmentranges = new AlignmentRangePB[num_of_chunks];

            bar = new shared(ProgressBar)();

            foreach (i; 0 .. num_of_chunks) {
                weights[i] = std.file.getSize(tmpfiles[i]); // set file size as weight
            }

            normalize(weights);

            foreach (i, ref range; alignmentranges) {
                auto bamfile = BamFile(tmpfiles[i]);
                bamfile.setBufferSize(memory_limit / 2 / num_of_chunks);
                range = bamfile.alignmentsWithProgress(
                // WTF is going on here? See this thread:
                // http://forum.dlang.org/thread/mailman.112.1341467786.31962.digitalmars-d@puremagic.com
                        (size_t j) { 
                            return (lazy float progress) { 
                                atomicStore(merging_progress[j], progress);
                                synchronized (bar) {
                                    bar.update(dotProduct(merging_progress, weights));
                                }
                            };
                        }(i));
            }

            writeBAM(stream, header_text, bam.reference_sequences,
                     nWayUnion!compareAlignmentCoordinates(alignmentranges),
                     compression_level,
                     task_pool);

            bar.finish();
        } else {
            alias ReturnType!(BamFile.alignments!withoutOffsets) AlignmentRange;
            auto alignmentranges = new AlignmentRange[num_of_chunks];

            foreach (i, ref range; alignmentranges) {
                auto bamfile = BamFile(tmpfiles[i]);
                bamfile.setBufferSize(memory_limit / 2 / num_of_chunks);
                range = bamfile.alignments!withoutOffsets;
            }

            writeBAM(stream, header_text, bam.reference_sequences,
                     nWayUnion!compareAlignmentCoordinates(alignmentranges),
                     compression_level,
                     task_pool);
        }

        return 0;
    } catch (Throwable e) {
        stderr.writeln("sambamba-sort: ", e.msg);
        return 1;
    }

    return 0;
}

/// parses \d+K, \d+M, \d+G
size_t parseMemory(string str) {
    auto initial_str = str.idup;
    size_t sz = 0;
    while (!str.empty) {
        if (isDigit(str[0]))
            sz *= 10, sz += str[0] - '0';
        else
            break;
        str = str[1 .. $];
    }
    if (str.empty) return sz;
    switch (str) {
        case "K":
        case "KiB":
            return sz << 10;
        case "KB":
            return sz * 1_000;
        case "M":
        case "MiB":
            return sz << 20;
        case "MB":
            return sz * 1_000_000;
        case "G":
        case "GiB":
            return sz << 30;
        case "GB":
            return sz * 1_000_000_000;
        default:
            throw new Exception("couldn't parse ", initial_str);
    }
}

/// Base name of the file corresponding to chunk number $(D chunk_num)
///
/// Params:
///     unsorted_fn - filename of unsorted BAM
///     chunk_num   - 0-based index of the chunk
///                                               
string chunkBaseName(string unsorted_fn, size_t chunk_num) {
    return baseName(unsorted_fn) ~ "." ~ to!string(chunk_num);
}

auto chunks(R)(R reads, size_t size_in_bytes) {

    static size_t approxSize(Alignment read) {
        return Alignment.sizeof * 3  / 2 + read.size_in_bytes;
    }

    // false means that each chunk may contain reads with different ref. ID
    return AlignmentRangeSplitter!(R, approxSize)(reads, size_in_bytes, false); 
}
