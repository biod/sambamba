/*
    This file is part of BioD.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
/// Writing a script/tool for processing BAM data often starts this way:
///
/// ------------------------
/// import bio.std.hts.bam.reader;
///
/// void main(string[] args) {
///     auto bam = new BamReader(args[1]); // open BAM file
///     foreach (read; bam.reads) {        // iterate through its reads
///         if (read.is_unmapped)
///             continue;                  // maybe skip unmapped ones
///         ...
///     }
/// }
/// ------------------------
///
/// Or, if a specific interval on the reference sequence is to be explored:
/// ------------------------
/// import bio.std.hts.bam.pileup;
/// ...
/// auto reads = bam["chr7"][50_000 .. 60_000]; // BAI index is required
/// foreach (column; makePileup(reads)) { ... } // see $(PMODULE pileup) docs
/// ------------------------
module bio.std.hts.bam.reader;

import bio.std.hts.bam.abstractreader;
public import bio.std.hts.sam.header;
public import bio.std.hts.bam.reference;
public import bio.std.hts.bam.region;
public import bio.std.hts.bam.read;
public import bio.std.hts.bam.tagvalue;
public import bio.std.hts.bam.readrange;

import bio.std.hts.bam.randomaccessmanager;
import bio.std.hts.bam.baifile;
import bio.std.hts.bam.bai.indexing;

import bio.core.utils.range;
import bio.core.utils.stream;
import bio.core.bgzf.inputstream;
public import bio.core.bgzf.virtualoffset;

import std.system;
import std.stdio;
import std.algorithm;
import std.range;
import std.conv;
import std.exception;
import std.parallelism;
import std.array;
import core.stdc.stdio;
import std.string;
import std.stdio;

/**
  BAM file reader, featuring parallel decompression of BGZF blocks.
 */
class BamReader : IBamSamReader {

    /**
      Creates reader associated with file or stream.
      (If stream constructor is used, no random access is possible.)
      $(BR)
      Optionally, task pool can be specified.
      It will be used to unpack BGZF blocks in parallel.

      Example:
      -------------------------------------------
      import std.parallelism, bio.std.hts.bam.reader;
      void main() {
        auto pool = new TaskPool(4); // use 4 threads
        scope (exit) pool.finish();  // don't forget!
        auto bam = new BamReader("file.bam", pool);
        ...
      }
      -------------------------------------------
     */
    this(contrib.undead.stream.Stream stream,
         TaskPool task_pool = taskPool) {
        _source_stream = new EndianStream(stream, Endian.littleEndian);
        _task_pool = task_pool;

        if (stream.seekable) {
            _stream_is_seekable = true;
        }

        initializeStreams();

        auto magic = _bam.readString(4);

        enforce(magic == "BAM\1", "Invalid file format: expected BAM\\1");

        readSamHeader();
        readReferenceSequencesInfo();

        // right after construction, we are at the beginning
        //                           of the list of reads

        if (_stream_is_seekable) {
            _reads_start_voffset = _decompressed_stream.virtualTell();
        }
    }

    /// ditto
    this(string filename, std.parallelism.TaskPool task_pool) {

        _filename = filename;
        _source_stream = getNativeEndianSourceStream();
        this(_source_stream, task_pool);
    }


    /// ditto
    this(string filename) {
      this(filename, std.parallelism.taskPool);
    }

    /**
      True if BAI file was found for this BAM file.
      This is necessary for any random-access operations.
      $(BR)
      Looks for files in the same directory which filename
      is either the file name of BAM file with '.bai' appended,
      or with the last extension replaced with '.bai'
      (that is, for $(I file.bam) paths $(I file.bai) and
      $(I file.bam.bai) will be checked)
     */
    bool has_index() @property {
        return _random_access_manager.found_index_file;
    }

    /**
      Creates BAI file. If $(I overwrite) is false, it won't touch
      existing index if it is already found.
     */
    void createIndex(bool overwrite = false) {
        if (has_index && !overwrite)
            return;
        Stream stream = new BufferedFile(filename ~ ".bai", FileMode.OutNew);
        scope(exit) stream.close();
        bio.std.hts.bam.bai.indexing.createIndex(this, stream);
        _bai_status = BaiStatus.notInitialized;
        _rndaccssmgr = null;
    }

    /** Filename, if the object was created via file name constructor,
        $(D null) otherwise.
     */
    string filename() @property const {
        return _filename;
    }

    /// If file ends with EOF block, returns virtual offset of the start of EOF block.
    /// Otherwise, returns virtual offset of the physical end of file.
    bio.core.bgzf.virtualoffset.VirtualOffset eofVirtualOffset() {
        return _random_access_manager.eofVirtualOffset();
    }

    /// Get BGZF block at a given file offset.
    bio.core.bgzf.block.BgzfBlock getBgzfBlockAt(ulong offset) {
        return _random_access_manager.getBgzfBlockAt(offset);
    }

    /**
      Returns: SAM header of the BAM file
     */
    bio.std.hts.sam.header.SamHeader header() @property {
        if (_header is null) {
            synchronized {
                if (_header is null) {
                    _header = new SamHeader(_headertext);
                    _headertext = null;
                }
            }
        }
        return _header;
    }

    /**
        Returns: information about reference sequences
     */
    const(bio.std.hts.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const nothrow {
        return _reference_sequences;
    }

    /**
        Range of all alignment records in the file.
        $(BR)
        Element type of the returned range depends on the policy.
        Default one is $(DPREF2 bam, readrange, withoutOffsets),
        in this case range element type is $(DPREF2 bam, read, BamRead).
        $(BR)
        The other option is $(DPREF2 bam, readrange, withOffsets),
        which allows to track read virtual offsets in the file.
        In this case range element type is $(DPREF2 bam, readrange, BamReadBlock).

        Example:
        ----------------------------------
        import bio.std.hts.bam.readrange;
        ...
        auto bam = new BamReader("file.bam");
        auto reads = bam.reads!withOffsets();
        writeln(reads.front.start_virtual_offset);
        ----------------------------------
     */
    auto reads(alias IteratePolicy=bio.std.hts.bam.readrange.withoutOffsets)() @property {
        auto _decompressed_stream = getDecompressedBamReadStream();
        return bamReadRange!IteratePolicy(_decompressed_stream, this);
    }

    static struct ReadsWithProgressResult(alias IteratePolicy, R, S) {
        this(R range, S stream,
             void delegate(lazy float p) progressBarFunc,
             void delegate() finishFunc)
        {
            _range = range;
            _stream = stream;
            _progress_bar_func = progressBarFunc;
            _finish_func = finishFunc;
        }

        static if (__traits(identifier, IteratePolicy) == "withOffsets") {
            auto front() @property {
                return _range.front;
            }
        } else static if (__traits(identifier, IteratePolicy) == "withoutOffsets") {
            auto front() @property {
                return _range.front.read;
            }
        } else static assert(0, __traits(identifier, IteratePolicy));

        bool empty() @property {
            auto result = _range.empty;
            if (_finish_func !is null && !_called_finish_func && result) {
                _called_finish_func = true;
                _finish_func();
            }
            return result;
        }

        void popFront() {
            _bytes_read += _range.front.read.size_in_bytes;
            _range.popFront();
            if (_progress_bar_func !is null) {
                _progress_bar_func(min(1.0,
                    cast(float)_bytes_read / (_stream.total_compressed_size *
                                              _stream.average_compression_ratio)));
            }
        }

        private R _range;
        private S _stream;
        private size_t _bytes_read;
        private void delegate(lazy float p) _progress_bar_func;
        private void delegate() _finish_func;
        private bool _called_finish_func = false;
    }

    /**
        Returns: range of all reads in the file, calling $(I progressBarFunc)
                 for each read.
        $(BR)
        $(I progressBarFunc) will be called
        each time next alignment is read, with the argument being a number from [0.0, 1.0],
        which is estimated progress percentage.
        $(BR)
        Notice that $(I progressBarFunc) takes $(D lazy) argument,
        so that the number of relatively expensive float division operations
        can be controlled by user.

        Once the iteration is finished (call to $(D empty) returned true),
        $(I finishFunc) will be called if provided.

        Example:
        ------------------------------------
        import std.functional, std.stdio, bio.std.hts.bam.reader;
        void progress(lazy float p) {
            static uint n;
            if (++n % 63 == 0) writeln(p); // prints progress after every 63 records
        }
        ...
        foreach (read; bam.readsWithProgress(toDelegate(&progress))) {
            ...
        }
        ------------------------------------
    */
    auto readsWithProgress(alias IteratePolicy=bio.std.hts.bam.readrange.withoutOffsets)
        (void delegate(lazy float p) progressBarFunc,
         void delegate() finishFunc=null)
    {
        auto _decompressed_stream = getDecompressedBamReadStream();
        auto reads_with_offsets = bamReadRange!withOffsets(_decompressed_stream, this);

        alias ReadsWithProgressResult!(IteratePolicy,
                       typeof(reads_with_offsets), BgzfInputStream) Result;

        return Result(reads_with_offsets, _decompressed_stream,
                      progressBarFunc, finishFunc);
    }

    ///
    void assumeSequentialProcessing() {
        _seqprocmode = true;
    }

    /// Part of IBamSamReader interface
    std.range.InputRange!(bio.std.hts.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads!withoutOffsets());
    }

    /**
      Returns: the read which starts at a given virtual offset.
     */
    bio.std.hts.bam.read.BamRead getReadAt(bio.core.bgzf.virtualoffset.VirtualOffset offset) {
        enforce(_random_access_manager !is null);
        return _random_access_manager.getReadAt(offset);
    }

    /**
      Returns: all reads located between two virtual offsets in the BAM file.

      $(BR)
      First offset must point to the start of an alignment record,
      and be strictly less than the second one.
      $(BR)
      For decompression, the task pool specified at the construction is used.
     */
    auto getReadsBetween(bio.core.bgzf.virtualoffset.VirtualOffset from,
                         bio.core.bgzf.virtualoffset.VirtualOffset to) {
        enforce(from <= to, "First offset must be less than second");
        enforce(_stream_is_seekable, "Stream is not seekable");

        return _random_access_manager.getReadsBetween(from, to);
    }

    /**
      Returns: all reads overlapping any region from a set.
    */
    auto getReadsOverlapping(BamRegion[] regions) {
	return _random_access_manager.getReads(regions);
    }

    /**
      Unmapped reads, i.e. reads at the end of file whose reference id is -1.
      The file must be coordinate-sorted and indexed.
     */
    auto unmappedReads() {
        enforce(_random_access_manager !is null);
        auto bai = _random_access_manager.getBai();

        VirtualOffset start;
        start = eofVirtualOffset();

        auto all_reads = this.reads();
        if (!all_reads.empty && all_reads.front.ref_id == -1)
            start = _reads_start_voffset;

        auto ioffsets = bai.indices[0 .. reference_sequences.length].retro()
                           .map!(index => index.ioffsets.retro()).joiner();
        if (!ioffsets.empty)
            start = ioffsets.front;

        auto stream = _random_access_manager.createStreamStartingFrom(start);
        auto r = bamReadRange!withOffsets(stream, this);
        while (!r.empty && r.front.ref_id != -1)
            r.popFront();
        return r;
    }

    /**
      Get BAI chunks containing all reads that overlap specified region.
      For $(I ref_id) = -1, use $(D unmappedReads) method.
     */
    bio.core.bgzf.chunk.Chunk[] getChunks(uint ref_id, int beg, int end) {
        enforce(_random_access_manager !is null);
        enforce(beg < end);

        return _random_access_manager.getChunks(BamRegion(ref_id, beg, end));
    }

    /**
      Returns reference sequence with id $(I ref_id).
     */
    bio.std.hts.bam.reference.ReferenceSequence reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return ReferenceSequence(_random_access_manager,
                                 ref_id,
                                 _reference_sequences[ref_id]);
    }

    /**
      Returns reference sequence named $(I ref_name).

      Example:
      ---------------------------
      import std.stdio, bio.std.hts.bam.reader;
      ...
      auto bam = new BamReader("file.bam");
      writeln(bam["chr2"].length);
      ---------------------------
     */
    bio.std.hts.bam.reference.ReferenceSequence opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " does not exist");
        auto ref_id = _reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    /**
      Check if reference named $(I ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        return null != (ref_name in _reference_sequence_dict);
    }

    /**
      Set buffer size for I/O operations. Values less than 4096 are disallowed.
      $(BR)
      This can help in multithreaded applications when several files are read
      simultaneously (e.g. merging).
     */
    void setBufferSize(size_t buffer_size) {
        enforce(buffer_size >= 4096, "Buffer size must be >= 4096 bytes");
        _buffer_size = buffer_size;
        _random_access_manager.setBufferSize(buffer_size);
    }

    package bool _seqprocmode; // available for bio.std.hts.bam.readrange;

private:

    string _filename;                       // filename (if available)
    Stream _source_stream;                  // compressed
    BgzfInputStream _decompressed_stream;   // decompressed
    Stream _bam;                            // decompressed + endian conversion
    bool _stream_is_seekable;

    // Virtual offset at which alignment records start.
    VirtualOffset _reads_start_voffset;

    BaiFile _dont_access_me_directly_use_bai_file_for_that;
    enum BaiStatus {
        notInitialized,
        initialized,
        fileNotFound
    }
    BaiStatus _bai_status = BaiStatus.notInitialized;

    void initBai() {
        if (_bai_status == BaiStatus.notInitialized) {
            synchronized {
                try {
                    _dont_access_me_directly_use_bai_file_for_that = BaiFile(_filename);
                    _bai_status = BaiStatus.initialized;
                } catch (Exception e) {
                    _bai_status = BaiStatus.fileNotFound;
                }
            }
        }
    }

    // provides access to index file
    @property ref BaiFile _bai_file() { // initialized lazily
        initBai();
        return _dont_access_me_directly_use_bai_file_for_that;
    };

    RandomAccessManager _rndaccssmgr; // unreadable for a purpose
    @property RandomAccessManager _random_access_manager() {
        if (_rndaccssmgr is null) {
            synchronized {
                initBai();

                if (_bai_status == BaiStatus.initialized) {
                    _rndaccssmgr = new RandomAccessManager(this, _bai_file);
                } else {
                    _rndaccssmgr = new RandomAccessManager(this);
                }

                _rndaccssmgr.setTaskPool(_task_pool);
                _rndaccssmgr.setBufferSize(_buffer_size);
            }
        }
        return _rndaccssmgr;
    }

    SamHeader _header;
    string _headertext; // for lazy SAM header parsing
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict; /// name -> index mapping

    TaskPool _task_pool;
    size_t _buffer_size = 4096; // buffer size to be used for I/O

    Stream getNativeEndianSourceStream() {
        assert(_filename !is null);
        Stream file = new bio.core.utils.stream.File(_filename);
        return new BufferedStream(file, _buffer_size);
    }

    Stream getSeekableCompressedStream() {
        if (_stream_is_seekable) {
            if (_filename !is null) {
                auto file = getNativeEndianSourceStream();
                version(development)
                {
                    std.stdio.stderr.writeln("[info] file size: ", file.size);
                }
                return new EndianStream(file, Endian.littleEndian);
            } else {
                _source_stream.seekSet(0);
                return _source_stream;
            }
        } else {
            return null;
        }
    }

    // get decompressed stream out of compressed BAM file
    BgzfInputStream getDecompressedStream() {
        auto compressed_stream = getSeekableCompressedStream();

        auto block_supplier = new StreamSupplier(compressed_stream is null ?
                                                 _source_stream :
                                                 compressed_stream);

        return new BgzfInputStream(block_supplier, _task_pool,
                                   _buffer_size);
    }

    // get decompressed stream starting from the first alignment record
    BgzfInputStream getDecompressedBamReadStream() {
        auto compressed_stream = getSeekableCompressedStream();

        if (compressed_stream !is null) {
            enforce(_reads_start_voffset != 0UL);

            compressed_stream.seekCur(_reads_start_voffset.coffset);
            auto block_supplier = new StreamSupplier(compressed_stream);
            auto stream = new BgzfInputStream(block_supplier, _task_pool,
                                              _buffer_size);
            stream.readString(_reads_start_voffset.uoffset);
            return stream;
        } else {
            // must be initialized in initializeStreams()
            return _decompressed_stream;
        }
    }

    // sets up the streams and ranges
    void initializeStreams() {

        _decompressed_stream = getDecompressedStream();
        _bam = new EndianStream(_decompressed_stream, Endian.littleEndian);
    }

    // initializes _header
    void readSamHeader() {
        int l_text;
        _bam.read(l_text);

        _headertext = cast(string)(_bam.readString(l_text));
    }

    // initialize _reference_sequences
    void readReferenceSequencesInfo() {
        int n_ref;
        _bam.read(n_ref);
        _reference_sequences = new ReferenceSequenceInfo[n_ref];
        foreach (i; 0..n_ref) {
            _reference_sequences[i] = ReferenceSequenceInfo(_bam);

            // provide mapping Name -> Index
            _reference_sequence_dict[_reference_sequences[i].name] = i;
        }
    }
}
