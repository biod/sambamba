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
import undead.stream;
import std.stdio;
import std.typecons;
import core.atomic;

import core.stdc.stdlib;
import core.stdc.string;

import sambamba.utils.common.progressbar;
import sambamba.utils.common.overwrite;
import sambamba.utils.common.tmpdir;
import sambamba.utils.common.filtering;
import sambamba.utils.common.file;

import thirdparty.mergesort;

void printUsage() {
    stderr.writeln("Usage: sambamba-sort [options] <input.bam>");
    stderr.writeln();
    stderr.writeln("Options: -m, --memory-limit=LIMIT");
    stderr.writeln("               approximate total memory limit for all threads (by default 2GB)");
    stderr.writeln("         --tmpdir=TMPDIR");
    stderr.writeln("               directory for storing intermediate files; default is system directory for temporary files");
    stderr.writeln("         -o, --out=OUTPUTFILE");
    stderr.writeln("               output file name; if not provided, the result is written to a file with .sorted.bam extension");
    stderr.writeln("         -n, --sort-by-name");
    stderr.writeln("               sort by read name instead of coordinate (lexicographical order)");
    stderr.writeln("         --sort-picard");
    stderr.writeln("               sort by query name like in picard");
    stderr.writeln("         -N, --natural-sort");
    stderr.writeln("               sort by read name instead of coordinate (so-called 'natural' sort as in samtools)");
    stderr.writeln("         -l, --compression-level=COMPRESSION_LEVEL");
    stderr.writeln("               level of compression for sorted BAM, from 0 to 9");
    stderr.writeln("         -u, --uncompressed-chunks");
    stderr.writeln("               write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("               show progressbar in STDERR");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("               use specified number of threads");
    stderr.writeln("         -F, --filter=FILTER");
    stderr.writeln("               keep only reads that satisfy FILTER");
}

version(standalone) {
    int main(string[] args) {
        return sort_main(args);
    }
}

private __gshared bool sort_by_name;
private __gshared bool natural_sort;
private __gshared bool picard_sort;
private bool show_progress;

private shared(ProgressBar) bar;
private shared(float[]) weights;
private shared(float[]) merging_progress;

class Sorter {

    BamReader bam;
    TaskPool task_pool;
    size_t memory_limit = 2u * 1024 * 1024 * 1024;
    string tmpdir;
    int compression_level = -1;
    bool uncompressed_chunks = false;
    string output_filename = null;
    string filename = null;
    string filter_str = null;

    this() {
        tmpdir = defaultTmpDir();
    }

    struct UnsortedChunk {
        size_t max_sz;

        this(size_t max_total_size) {
            max_sz = max_total_size;
            while (read_storage is null && max_sz > 65536) {
                read_storage = cast(ubyte*)core.stdc.stdlib.malloc(max_sz);
                if (read_storage is null)
                  max_sz /= 2;
            }
            _reads_capa = 1024;
            auto sz = BamRead.sizeof * _reads_capa;
            _reads = cast(BamRead*)core.stdc.stdlib.malloc(sz);
            if (_reads is null) {
                throw new Exception("alloc failed: no space for read pointers");
            }
        }

        void clear() {
            _used = 0;
            _n_reads = 0;
            if (_low_memory) {
                auto realloc_storage = cast(ubyte*)core.stdc.stdlib.realloc(read_storage, max_sz / 2);
                if (realloc_storage !is null) {
                    max_sz /= 2;
                    read_storage = realloc_storage;
                    stderr.writeln("reduced maximum buffer size to ", max_sz);
                }
                _low_memory = false;
            }
        }

        void fill(R)(R* reads) {
            while (!reads.empty) {
                auto read = reads.front;
                auto len = read.raw_data.length;
                if (len + _used > max_sz)
                    break;

                if (_n_reads == _reads_capa) {
                    auto realloc_reads = cast(BamRead*)core.stdc.stdlib.realloc(_reads, 2 * _reads_capa * BamRead.sizeof);
                    if (realloc_reads is null) {
                        _low_memory = true;
                        stderr.writeln("realloc failed: system low on memory, limited to ", _reads_capa, " reads in buffer");
                        break;
                    } else {
                        _reads_capa *= 2;
                        _reads = realloc_reads;
                    }
                }
                core.stdc.string.memcpy(read_storage + _used, read.raw_data.ptr, len);
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
                auto realloc_storage = cast(ubyte*)core.stdc.stdlib.realloc(read_storage, len);
                if (realloc_storage is null) {
                    throw new Exception("realloc failed: not enough memory for read");
                } else {
                    read_storage = realloc_storage;
                    max_sz = len;
                }
                _used = len;
                read_storage[0 .. len] = read.raw_data[];
                _reads[0].raw_data = read_storage[0 .. _used];
                _reads[0].associateWithReader(read.reader);
                reads.popFront();
            }
        }

        void free() {
            core.stdc.stdlib.free(read_storage);
            core.stdc.stdlib.free(_reads);
        }

        BamRead[] reads() @property { return _reads[0 .. _n_reads]; }
        BamRead* _reads;
        size_t _reads_capa;
        size_t _n_reads;

        ubyte* read_storage;
        size_t _used;
        bool _low_memory;
    }

    static BamRead[] sortChunk(size_t n, BamRead[] chunk, TaskPool task_pool) {
        version (development) {
        StopWatch sw;
        sw.start();
        stderr.writeln("Sorting chunk #", n, "...");
        }
        auto buf = cast(BamRead*)core.stdc.stdlib.malloc(chunk.length * BamRead.sizeof);
        BamRead[] tmp = buf[0 .. chunk.length];
        scope (exit) core.stdc.stdlib.free(buf);
        if (sort_by_name) {
            mergeSort!(compareReadNames, false)(chunk, task_pool, tmp);
        } else if (natural_sort) {
            mergeSort!(mixedCompareReadNames, false)(chunk, task_pool, tmp);
        } else if (picard_sort) {
            mergeSort!(compareReadNamesAsPicard, false)(chunk, task_pool, tmp);
        } else {
            mergeSort!(compareCoordinatesAndStrand, false)(chunk, task_pool, tmp);
        }
        version (development) {
        stderr.writeln("Finished sorting of chunk #", n, " in ", sw.peek().seconds, "s");
        }
        return chunk;
    }

    void sort() {
        auto filter = createFilterFromQuery(filter_str);

        createHeader();

        size_t buf_size = 16_000_000;
        bam.setBufferSize(buf_size);
        bam.assumeSequentialProcessing();

        auto output_buffer_ptr = cast(ubyte*)core.stdc.stdlib.malloc(buf_size);
        output_buffer = output_buffer_ptr[0 .. buf_size];
        scope(exit)core.stdc.stdlib.free(output_buffer_ptr);

        if (show_progress) {
            stderr.writeln("Writing sorted chunks to temporary directory...");
            bar = new shared(ProgressBar)();
            auto reads = bam.readsWithProgress(
                                  (lazy float p){ bar.update(p); }
                              );
            auto filtered_reads = filtered(reads, filter);
            writeSortedChunks(filtered_reads);
            bar.finish();
        } else {
            auto filtered_reads = filtered(bam.reads!withoutOffsets(), filter);
            writeSortedChunks(filtered_reads);
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
        else if (natural_sort)
            mergeSortedChunks!mixedCompareReadNames();
        else if (picard_sort)
            mergeSortedChunks!compareReadNamesAsPicard();
        else
            mergeSortedChunks!compareCoordinatesAndStrand();
    }

    private void createHeader() {
        header = bam.header;
        header.sorting_order = (sort_by_name || natural_sort || picard_sort) ? SortingOrder.queryname :
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

        version(development) stderr.writeln("waiting for the last sorting task...");

        if (sorting_task !is null)
            dump(sorting_task.yieldForce());

        // handle empty BAM file case
        if (k == 0) {
            ++k;
            dump([]);
        }
    }

    private ubyte[] output_buffer;

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

        if (k <= 1) {
            level = compression_level;
            fn = output_filename;
        } else {
            fn = tmpFile(chunkBaseName(filename, num_of_chunks), tmpdir);
            tmpfiles ~= fn;
        }

        auto stream = bufferedFile(fn, FileMode.OutNew, 0);
        stream.buffer = output_buffer;
        scope(failure) stream.close();

        // if there's only one chunk, we will write straight to the output file
        auto writer = scoped!BamWriter(stream, level, task_pool, 32_000_000);
        writer.setFilename(fn);

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

        if (num_of_chunks <= 1) {
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
        auto stream = bufferedFile(output_filename, FileMode.OutNew,
                                   output_buf_size);
        scope(failure) stream.close();

        alias ReturnType!(BamReader.readsWithProgress!withoutOffsets) AlignmentRangePB;
        auto alignmentranges = new AlignmentRangePB[num_of_chunks];

        if (show_progress) {
            stderr.writeln("Merging sorted chunks...");
            float[] weights1;
            weights1.length = num_of_chunks;
            merging_progress.length = num_of_chunks;
            merging_progress[] = 0.0;

            bar = new shared(ProgressBar)();

            foreach (i; 0 .. num_of_chunks) {
                weights1[i] = std.file.getSize(tmpfiles[i]); // set file size as weight
            }

            normalize(weights1);
            weights = cast(shared)weights1.dup();
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

        auto writer = scoped!BamWriter(stream, compression_level, task_pool,
                                       2 * output_buf_size);
        writer.setFilename(output_filename);
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

        auto sorter = scoped!Sorter();

        getopt(args,
               std.getopt.config.caseSensitive,
               "memory-limit|m",        &memory_limit_str,
               "tmpdir",                &sorter.tmpdir,
               "out|o",                 &sorter.output_filename,
               "sort-by-name|n",        &sort_by_name,
               "natural-sort|N",        &natural_sort,
               "sort-picard",           &picard_sort,
               "uncompressed-chunks|u", &sorter.uncompressed_chunks,
               "compression-level|l",   &sorter.compression_level,
               "show-progress|p",       &show_progress,
               "nthreads|t",            &n_threads,
               "filter|F",              &sorter.filter_str);

        if ((sort_by_name && (natural_sort || picard_sort)) || (natural_sort && picard_sort)) {
            stderr.writeln("only one of -n and -N and -s parameters can be provided");
            return -1;
        }

        if (sorter.output_filename is null) {
            sorter.output_filename = setExtension(args[1], "sorted.bam");
        }

        protectFromOverwrite(args[1], sorter.output_filename);
        sorter.tmpdir = randomSubdir(sorter.tmpdir);

        if (memory_limit_str !is null) {
            sorter.memory_limit = parseMemory(memory_limit_str);
            if (sorter.memory_limit / max(n_threads, 1) < 100_000) {
                throw new Exception("memory limit per thread can't be less than 100Kb");
            }
        }

        sorter.memory_limit = (sorter.memory_limit * 5) / 6;

        sorter.task_pool = new TaskPool(n_threads);
        scope(exit) sorter.task_pool.finish();
        sorter.bam = new BamReader(args[1], sorter.task_pool);
        sorter.filename = args[1];

        import core.memory;
        GC.disable();

        sorter.sort();

        try {
          std.file.rmdirRecurse(sorter.tmpdir);
        } catch (FileException e) {
          // Ignore errors removing temporary directories, due to NFS failure under load
          // https://github.com/chapmanb/bcbio-nextgen/issues/784
        }

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
