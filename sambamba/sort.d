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
import std.datetime;
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

import std.c.stdlib;
import std.c.string;

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

private __gshared bool sort_by_name;
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

    struct UnsortedChunk {
        size_t max_sz;

        this(size_t max_total_size) {
            max_sz = max_total_size;
            read_storage = cast(ubyte*)std.c.stdlib.malloc(max_sz);
            _reads_capa = 1024;
            auto sz = BamRead.sizeof * _reads_capa;
            _reads = cast(BamRead*)std.c.stdlib.malloc(sz);
        }

        void clear() {
            _used = 0;
            _n_reads = 0;
        }

        void fill(R)(R* reads) {
            while (!reads.empty) {
                auto read = reads.front;
                auto len = read.raw_data.length;
                if (len + _used > max_sz)
                    break;
                
                std.c.string.memcpy(read_storage + _used, read.raw_data.ptr, len);
                if (_n_reads == _reads_capa) {
                    _reads_capa *= 2;
                    _reads = cast(BamRead*)std.c.stdlib.realloc(_reads, _reads_capa * BamRead.sizeof);
                }
                _reads[_n_reads].raw_data = read_storage[_used .. _used + len];
                _reads[_n_reads].associateWithReader(read.reader);

                _n_reads += 1;
                _used += len;

                reads.popFront();
            }

            if (_n_reads == 0) {
                auto read = reads.front;
                auto len = read.raw_data.length;
                assert(len > max_sz);
                _n_reads = 1;
                read_storage = cast(ubyte*)std.c.stdlib.realloc(read_storage, len);
                _used = len;
                read_storage[0 .. len] = read.raw_data[];
                _reads[0].raw_data = read_storage[0 .. _used];
                _reads[0].associateWithReader(read.reader);
                reads.popFront();
            }
        }

        void free() {
            std.c.stdlib.free(read_storage);
            std.c.stdlib.free(_reads);
        }

        BamRead[] reads() @property { return _reads[0 .. _n_reads]; }
        BamRead* _reads;
        size_t _reads_capa;
        size_t _n_reads;

        ubyte* read_storage;
        size_t _used;
    }

    static BamRead[] sortChunk(size_t n, BamRead[] chunk, TaskPool task_pool) {
        version (development) {
        StopWatch sw;
        sw.start();
        stderr.writeln("Sorting chunk #", n, "...");
        }
        auto buf = cast(BamRead*)std.c.stdlib.malloc(chunk.length * BamRead.sizeof);
        BamRead[] tmp = buf[0 .. chunk.length];
        scope (exit) std.c.stdlib.free(buf);
        if (!sort_by_name) {
            mergeSort!(compareCoordinates, false)(chunk, task_pool, tmp);
        } else {
            mergeSort!(compareReadNames, false)(chunk, task_pool);
        }
        version (development) {
        stderr.writeln("Finished sorting of chunk #", n, " in ", sw.peek().seconds, "s");
        }
        return chunk;
    }

    void sort() {
        createHeader();

        bam.setBufferSize(16_000_000);
        bam.assumeSequentialProcessing();
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

        scope(success) {
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

    private size_t k; // number of sorting tasks submitted

    private void writeSortedChunks(R)(R reads_) {
        auto buf1 = UnsortedChunk(memory_limit / 2);
        auto buf2 = UnsortedChunk(memory_limit / 2);
        scope(exit) { buf1.free(); buf2.free(); }

        version (development) {
        StopWatch sw;
        sw.start();
        scope(exit) {
            sw.stop();
            stderr.writeln("Wrote ", num_of_chunks, " sorted chunk",
                           num_of_chunks == 1 ? "" : "s",
                           " in ", sw.peek().seconds, " seconds");
        }
        }

        Task!(sortChunk, size_t, BamRead[], TaskPool)* sorting_task;

        auto reads = reads_;

        // 1) fill buf1
        // 2) sort buf1 in parallel with filling buf2
        // 3) dump buf1
        // 4) sort buf2 in parallel with filling buf1
        // ...
        while (!reads.empty) {
            buf1.clear();
            version(development) {
            stderr.writeln("Reading chunk #", k + 1);
            StopWatch sw_inner;
            sw_inner.start();
            }
            buf1.fill(&reads);

            version(development)
            stderr.writeln("Finished reading of chunk #", k + 1,
                           " in ", sw_inner.peek().seconds, "s");

            BamRead[] sorted_reads;
            if (sorting_task !is null)
                sorted_reads = sorting_task.yieldForce();

            sorting_task = task!sortChunk(k + 1, buf1.reads, task_pool);
            task_pool.put(sorting_task);
            ++k;

            if (sorted_reads.length > 0)
                dump(sorted_reads);

            swap(buf1, buf2);
        }

        if (sorting_task !is null)
            dump(sorting_task.yieldForce());
    }

    // if there's more than one chunk, first call to dump will be when k == 2
    // because we first submit sorting task and only when dump previous chunk
    private void dump(BamRead[] sorted_reads) {
        version(development) {
            stderr.writeln("Dumping chunk #", num_of_chunks + 1, " to disk...");
            StopWatch sw;
            sw.start();
        }

        int level = uncompressed_chunks ? 0 : 1;

        string fn;

        if (k == 1) { 
            level = compression_level;
            fn = output_filename;
        } else {
            fn = tmpFile(chunkBaseName(filename, num_of_chunks), tmpdir);
            tmpfiles ~= fn;
        }

        Stream stream = new BufferedFile(fn, FileMode.OutNew, 16_000_000);
        scope(failure) stream.close();

        // if there's only one chunk, we will write straight to the output file
        auto writer = new BamWriter(stream, level, task_pool, 32_000_000);

        writer.writeSamHeader(header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);

        foreach (read; sorted_reads)
            writer.writeRecord(read);

        writer.finish();

        version(development)
        stderr.writeln("Finished dumping of chunk #", num_of_chunks + 1,
                       " in ", sw.peek().seconds, "s");

        num_of_chunks += 1;
    }

    private void mergeSortedChunks(alias comparator)() {

        if (num_of_chunks == 1) {
            // dump() wrote it to destination already
            return;
        }

        version(development) {
        StopWatch sw;
        sw.start();
        scope(exit) {
            sw.stop();
            stderr.writeln("Merging took ", sw.peek().seconds, " seconds");
        }
        }

        auto input_buf_size = min(16_000_000, memory_limit / 4 / num_of_chunks);
        auto output_buf_size = min(64_000_000, memory_limit / 6);
        Stream stream = new BufferedFile(output_filename, FileMode.OutNew, 
                                         output_buf_size);
        scope(failure) stream.close();

        alias ReturnType!(BamReader.readsWithProgress!withoutOffsets) AlignmentRangePB;
        auto alignmentranges = new AlignmentRangePB[num_of_chunks];

        if (show_progress) {
            stderr.writeln("Merging sorted chunks...");
            weights.length = num_of_chunks;
            merging_progress.length = num_of_chunks;
            merging_progress[] = 0.0;

            bar = new shared(ProgressBar)();

            foreach (i; 0 .. num_of_chunks) {
                weights[i] = std.file.getSize(tmpfiles[i]); // set file size as weight
            }

            normalize(cast()weights);
        }

        foreach (i; 0 .. num_of_chunks) {
            auto bamfile = new BamReader(tmpfiles[i], task_pool);
            bamfile.setBufferSize(input_buf_size);
            bamfile.assumeSequentialProcessing();
            if (show_progress)
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
            else
                alignmentranges[i] = bamfile.readsWithProgress(null);
        }

        auto writer = new BamWriter(stream, compression_level, task_pool,
                2 * output_buf_size);
        scope(exit) writer.finish();
        writer.writeSamHeader(header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);

        foreach (read; nWayUnion!comparator(alignmentranges))
            writer.writeRecord(read);

        if (show_progress)
            bar.finish();
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

        sorter.memory_limit = (sorter.memory_limit * 2) / 3;

        sorter.task_pool = new TaskPool(n_threads);
        scope(exit) sorter.task_pool.finish();
        sorter.bam = new BamReader(args[1], sorter.task_pool);
        sorter.filename = args[1];

        sorter.sort();
        return 0;
    } catch (Throwable e) {
        version (development) {
            throw e;
        } else {
            stderr.writeln("sambamba-sort: ", e.msg);
        }
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
