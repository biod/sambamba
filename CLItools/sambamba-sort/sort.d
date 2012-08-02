module sambamba.sort;

import bamfile;
import bamoutput;
import samheader;
import alignment;

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
    writeln("Usage: sambamba-sort [options] <input.bam>");
    writeln();
    writeln("Options: -m, --memory-limit=LIMIT");
    writeln("               approximate memory limit (it's not guaranteed that it won't be exceeded, because of garbage collection)");
    writeln("         --tmpdir=TMPDIR");
    writeln("               directory for storing intermediate files; default is system directory for temporary files");
    writeln("         -o, --out=OUTPUTFILE");
    writeln("               output file name; if not provided, the result is written to a file with .sorted.bam extension");
    writeln("         -n, --sort-by-name");
    writeln("               sort by read name instead of coordinate");
    writeln("         -l, --compression-level=COMPRESSION_LEVEL");
    writeln("               level of compression for sorted BAM, from 0 to 9");
    writeln("         -u, --uncompressed-chunks");
    writeln("               write sorted chunks as uncompressed BAM (default is writing with compression level 1), that might be faster in some cases but uses more disk space");
    writeln("         -p, --show-progress");
    writeln("               show progressbar in STDERR");
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
        string output_filename = setExtension(args[1], "sorted.bam");
        bool uncompressed_chunks;
        int compression_level = -1;

        getopt(args,
               std.getopt.config.caseSensitive,
               "memory-limit|m",        &memory_limit_str,
               "tmpdir",                &tmpdir,
               "out|o",                 &output_filename,
               "sort-by-name|n",        &sort_by_name,
               "uncompressed-chunks|u", &uncompressed_chunks,
               "compression-level|l",   &compression_level,
               "show-progress|p",       &show_progress);

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
        Stream stream = new BufferedFile(output_filename, FileMode.Out,
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

/// Absolute path of temporary file.
///
/// Params:
///     filename - base name
///     tmpdir   - temporary directory
///                                   
string tmpFile(string filename, string tmpdir) {
    if (tmpdir != null) {
        return buildPath(tmpdir, filename);
    }
version(Windows)
    return buildPath(std.process.getenv("TEMP"), filename);
else version(Posix)
    return "/tmp/" ~ filename;
}

// Constructs range of chunks where total size of alignments
// in a chunk does not exceed given amount of bytes.
struct Chunks(R) 
    if (isInputRange!R && is(Unqual!(ElementType!R) == Alignment))
{
    this(R range, size_t size) {
        _range = range;
        _size = size;
        _appender = appender!(Alignment[])();
        getNextChunk();
    }

    private {
        R _range;
        bool _empty;
        size_t _size;
        Appender!(Alignment[]) _appender;
    }

    bool empty() @property {
        return _empty;
    }

    Alignment[] front() @property {
        return _appender.data.dup;
    }

    void popFront() {
        _appender.clear();
        getNextChunk(); 
    }

    private void getNextChunk() {
        if (_range.empty) {
            _empty = true;
            return;
        } 

        auto first_read = _range.front;
        _range.popFront();

        size_t total_size = first_read.size_in_bytes;
        auto average_size_estimate = total_size + Alignment.sizeof * 3/2;

        _appender.reserve(_size / average_size_estimate);
        _appender.put(first_read);

        while (total_size <= _size && !_range.empty) {
            auto read = _range.front;
            total_size += Alignment.sizeof * 3 / 2 + // for mergesort
                          read.size_in_bytes;
            _appender.put(read);
            _range.popFront();
        }
        debug {
            import std.stdio;
            writeln(total_size);
        }
    }
}

auto chunks(R)(R alignments, size_t chunk_size) {
    return Chunks!R(alignments, chunk_size);
}
