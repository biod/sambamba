import bamfile;
import bamoutput;
import samheader;
import alignment;

import std.range;
import std.algorithm;
import std.traits;
import std.array;
import std.ascii;
import std.parallelism;
import std.getopt;
import std.path;
import std.file;
import std.stream;
import std.stdio;

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
}

int main(string[] args) {

    if (args.length < 2) {
        printUsage();
        return 0;
    }

    try {
        string memory_limit_str = null;

        size_t memory_limit = 512 * 1024 * 1024;
        string tmpdir = null;
        string output_filename = setExtension(args[1], "sorted.bam");

        getopt(args,
               std.getopt.config.caseSensitive,
               "memory-limit|m", &memory_limit_str,
               "tmpdir",         &tmpdir,
               "out|o",          &output_filename);

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
        header.sorting_order = SortingOrder.coordinate;
        auto header_text = toSam(header);

        // ----------------------- sort chunks -------------------------------------

        static Alignment[] sortChunk(Alignment[] chunk) {
            mergeSort!compareAlignmentCoordinates(chunk, true); // threaded
            return chunk;
        }

        string[] tmpfiles;

        auto num_of_chunks = 0;
        foreach (chunk; map!sortChunk(chunks(bam.alignments, memory_limit)))
        {
            auto fn = tmpFile(chunkBaseName(args[1], num_of_chunks), tmpdir);
            tmpfiles ~= fn;

            Stream stream = new BufferedFile(fn, FileMode.Out);
            scope(exit) stream.close();

            writeBAM(stream, header_text, bam.reference_sequences, 
                     chunk, 1, task_pool);

            num_of_chunks += 1;
        }

        // ---------------------- merge sorted chunks ------------------------------

        alias ReturnType!(BamFile.alignments) AlignmentRange;
        auto alignmentranges = new AlignmentRange[num_of_chunks];

        // half of memory is for input buffers
        foreach (i, ref range; alignmentranges) {
            auto bamfile = BamFile(tmpfiles[i]);
            bamfile.setBufferSize(memory_limit / 2 / num_of_chunks);
            range = bamfile.alignments;
        }

        // and another half is for output buffers
        Stream stream = new BufferedFile(output_filename, FileMode.Out,
                                         memory_limit / 2);
        scope(exit) stream.close();

        writeBAM(stream, header_text, bam.reference_sequences,
                 nWayUnion!compareAlignmentCoordinates(alignmentranges),
                 -1,
                 task_pool);

        // ---------------- remove temporary files -----------------------------

        foreach (tmpfile; tmpfiles) {
            remove(tmpfile);
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
            return sz * 1024UL;
        case "KB":
            return sz * 1000UL;
        case "M":
        case "MiB":
            return sz * 1024UL * 1024UL;
        case "MB":
            return sz * 1000UL * 1000UL;
        case "G":
        case "GiB":
            return sz * 1024UL * 1024UL * 1024UL;
        case "GB":
            return sz * 1000UL * 1000UL * 1000UL;
        default:
            throw new Exception("couldn't parse ", initial_str);
    }
}

/// Comparison function for alignments
bool compareAlignmentCoordinates(Alignment a1, Alignment a2) {
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    if (a1.position < a2.position) return true;
    return false;
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

// copy-pasted from std.algorithm. The only difference is that front is not a ref.
import std.container;
import std.functional;

struct NWayUnion(alias less, RangeOfRanges)
{
    private alias .ElementType!(.ElementType!RangeOfRanges) ElementType;
    private alias binaryFun!less comp;
    private RangeOfRanges _ror;
    static bool compFront(.ElementType!RangeOfRanges a,
            .ElementType!RangeOfRanges b)
    {
        return comp(b.front, a.front);
    }
    BinaryHeap!(RangeOfRanges, compFront) _heap;

    this(RangeOfRanges ror)
    {
        _ror = remove!("a.empty", SwapStrategy.unstable)(ror);
        _heap.acquire(_ror);
    }

    @property bool empty() { return _ror.empty; }

    @property ElementType front() // <-------- the only difference
    {
        return _heap.front.front;
    }

    void popFront()
    {
        _heap.removeFront();
        _ror.back.popFront();
        if (_ror.back.empty)
        {
            _ror.popBack();
            return;
        }
        _heap.conditionalInsert(_ror.back) || assert(false);
    }
}

NWayUnion!(less, RangeOfRanges) nWayUnion
(alias less = "a < b", RangeOfRanges)
(RangeOfRanges ror)
{
    return typeof(return)(ror);
}
