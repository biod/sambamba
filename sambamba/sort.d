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
module sambamba.sort;

import bio.bam.reader;
import bio.bam.writer;
import bio.sam.header;
import bio.bam.read;
import bio.bam.splitter;
import bio.core.utils.tmpfile;

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

import sambamba.utils.common.progressbar;

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
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("               use specified number of threads");
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

class Sorter {

    BamReader bam;
    TaskPool task_pool;
    size_t memory_limit = 512 * 1024 * 1024;
    string tmpdir = null;
    int compression_level = -1;
    bool uncompressed_chunks = false;
    string output_filename = null;
    string filename = null;

    private static auto chunks(R)(R reads, size_t size_in_bytes) {

        static size_t approxSize(BamRead read) {
            return BamRead.sizeof * 3  / 2 + read.size_in_bytes;
        }

        // false means that each chunk may contain reads with different ref. ID
        return ReadRangeSplitter!(R, approxSize)(reads, size_in_bytes, false); 
    }

    static BamRead[] sortChunk(BamRead[] chunk, TaskPool task_pool) {
        if (!sort_by_name) {
            mergeSort!compareCoordinates(chunk, task_pool);
        } else {
            mergeSort!compareReadNames(chunk, task_pool);
        }
        return chunk;
    }

    void sort() {
        createHeader();

        if (show_progress) {
            stderr.writeln("Writing sorted chunks to temporary directory...");
            bar = new shared(ProgressBar)();
            auto reads = bam.readsWithProgress(
                                  (lazy float p){ bar.update(p); }
                              );
            writeSortedChunks(reads);
            bar.finish();
        } else {
            writeSortedChunks(bam.reads!withoutOffsets);
        }

        scope(exit) {
            if (tmpfiles.length > 1) {
                // if length == 1, the file was moved
                foreach (tmpfile; tmpfiles) {
                    remove(tmpfile);
                }
            }
        }

        if (sort_by_name)
            mergeSortedChunks!compareReadNames();
        else
            mergeSortedChunks!compareCoordinates();
    }

    private void createHeader() {
        header = bam.header;
        header.sorting_order = sort_by_name ? SortingOrder.queryname :
                                              SortingOrder.coordinate;
    }

    private void writeSortedChunks(R)(R reads) {
        auto unsorted_chunks = chunks(reads, memory_limit);
        foreach (unsorted_chunk; unsorted_chunks)
        {
            auto chunk = sortChunk(unsorted_chunk, task_pool);

            auto fn = tmpFile(chunkBaseName(filename, num_of_chunks), tmpdir);
            tmpfiles ~= fn;

            Stream stream = new BufferedFile(fn, FileMode.Out);
            scope(failure) stream.close();

            auto writer = new BamWriter(stream, 
                                        uncompressed_chunks ? 0 : 1,
                                        task_pool);
            scope(exit) writer.finish();

            writer.writeSamHeader(header);
            writer.writeReferenceSequenceInfo(bam.reference_sequences);

            foreach (read; chunk)
                writer.writeRecord(read);

            num_of_chunks += 1;
        }
    }

    private void mergeSortedChunks(alias comparator)() {

        // optimization: if there's only one chunk, just rename it
        if (num_of_chunks == 1) {
            std.file.rename(tmpfiles[0], output_filename);
            return;
        }

        // half of memory is for input buffers
        // and another half is for output buffers
        Stream stream = new BufferedFile(output_filename, FileMode.OutNew, 
                                         memory_limit / 2);
        scope(failure) stream.close();

        if (show_progress) {
            stderr.writeln("Merging sorted chunks...");
            weights.length = num_of_chunks;
            merging_progress.length = num_of_chunks;
            merging_progress[] = 0.0;

            alias ReturnType!(BamReader.readsWithProgress!withoutOffsets) AlignmentRangePB;
            auto alignmentranges = new AlignmentRangePB[num_of_chunks];

            bar = new shared(ProgressBar)();
            scope(exit) bar.finish();

            foreach (i; 0 .. num_of_chunks) {
                weights[i] = std.file.getSize(tmpfiles[i]); // set file size as weight
            }

            normalize(cast()weights);

            foreach (i; 0 .. num_of_chunks) {
                auto bamfile = new BamReader(tmpfiles[i], task_pool);
                bamfile.setBufferSize(memory_limit / 2 / num_of_chunks);
                alignmentranges[i] = bamfile.readsWithProgress(
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

            auto writer = new BamWriter(stream, compression_level, task_pool);
            scope(exit) writer.finish();
            writer.writeSamHeader(header);
            writer.writeReferenceSequenceInfo(bam.reference_sequences);

            foreach (read; nWayUnion!comparator(alignmentranges))
                writer.writeRecord(read);

        } else {
            alias ReturnType!(BamReader.reads!withoutOffsets) AlignmentRange;
            auto alignmentranges = new AlignmentRange[num_of_chunks];

            foreach (i; 0 .. num_of_chunks) {
                auto bamfile = new BamReader(tmpfiles[i]);
                bamfile.setBufferSize(memory_limit / 2 / num_of_chunks);
                alignmentranges[i] = bamfile.reads!withoutOffsets;
            }

            auto writer = new BamWriter(stream, compression_level, task_pool);
            scope(exit) writer.finish();
            writer.writeSamHeader(header);
            writer.writeReferenceSequenceInfo(bam.reference_sequences);

            foreach (read; nWayUnion!comparator(alignmentranges))
                writer.writeRecord(read);
        }
    }

    private {
        SamHeader header;

        string[] tmpfiles; // temporary file names
        size_t num_of_chunks; // number of temporary files
    }
}

int sort_main(string[] args) {

    if (args.length < 2) {
        printUsage();
        return 0;
    }

    try {
        string memory_limit_str = null;
        uint n_threads = totalCPUs;

        auto sorter = new Sorter();

        getopt(args,
               std.getopt.config.caseSensitive,
               "memory-limit|m",        &memory_limit_str,
               "tmpdir",                &sorter.tmpdir,
               "out|o",                 &sorter.output_filename,
               "sort-by-name|n",        &sort_by_name,
               "uncompressed-chunks|u", &sorter.uncompressed_chunks,
               "compression-level|l",   &sorter.compression_level,
               "show-progress|p",       &show_progress,
               "nthreads|t",            &n_threads);

        if (sorter.output_filename is null) {
            sorter.output_filename = setExtension(args[1], "sorted.bam");
        }

        if (memory_limit_str !is null) {
            sorter.memory_limit = parseMemory(memory_limit_str);
        }

        // TODO: find a better formula
        sorter.memory_limit /= 4;

        sorter.task_pool = new TaskPool(n_threads);
        scope(exit) sorter.task_pool.finish();
        sorter.bam = new BamReader(args[1], sorter.task_pool);
        sorter.filename = args[1];

        sorter.sort();
        return 0;
    } catch (Throwable e) {
        stderr.writeln("sambamba-sort: ", e.msg);
        return 1;
    }
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
