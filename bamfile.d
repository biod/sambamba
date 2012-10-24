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
module bamfile;

public import samheader;
public import reference;
public import alignment;
public import virtualoffset;
public import tagvalue;
public import alignmentrange;
import bgzfrange;
import chunkinputstream;
import randomaccessmanager;
import bai.read;
import bai.chunk;
import utils.range;

import utils.stream;
import std.system;
import std.stdio;
import std.algorithm : map, min;
import std.range : zip;
import std.conv : to;
import std.exception : enforce;
import std.parallelism;
import std.array : uninitializedArray;
import core.stdc.stdio;
import std.string;

/**
  Represents BAM file
 */
struct BamFile {

    this(Stream stream, TaskPool task_pool = taskPool) {
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
        //                           of the list of alignments

        if (_stream_is_seekable) {
            _alignments_start_voffset = _decompressed_stream.virtualTell();
        }
    }

    /**
      Constructor taking filename of BAM file to open,
      and optionally, task pool to use.
      
      Currently, opens the file read-only since library
      has no support for writing yet.
     */
    this(string filename, TaskPool task_pool = taskPool) {

        _filename = filename;
        _source_stream = getNativeEndianSourceStream();
        this(_source_stream, task_pool);
    }
  
    /// True if associated BAI file was found
    bool has_index() @property {
        return _random_access_manager.found_index_file;
    }

    /// If file ends with EOF block, returns virtual offset of the start of EOF block.
    /// Otherwise, returns virtual offset of the physical end of file.
    VirtualOffset eofVirtualOffset() {
        return _random_access_manager.eofVirtualOffset();
    }

    /// Get BGZF block at a given file offset.
    BgzfBlock getBgzfBlockAt(ulong offset) {
        return _random_access_manager.getBgzfBlockAt(offset);
    }

    /*
       Get SAM header of file.
     */
    SamHeader header() @property {
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
    ReferenceSequenceInfo[] reference_sequences() @property {
        return _reference_sequences;
    }

    /**
        Returns: range of all alignments in the file.

        However, using several ranges is not recommended since it can hurt
        disk access performance.
     */
    auto alignments(alias IteratePolicy=withoutOffsets)() @property {
        auto _decompressed_stream = getDecompressedAlignmentStream();
        return alignmentRange!IteratePolicy(_decompressed_stream);
    }

    /**
        Returns: range of all alignments in the file, calling $(D progressBarFunc)
                 for each alignment. 

        $(D progressBarFunc) will be called
        each time next alignment is read, with the argument being a number from [0.0, 1.0],
        which is estimated progress percentage.
    */
    auto alignmentsWithProgress(alias IteratePolicy=withoutOffsets)
        (void delegate(lazy float p) progressBarFunc) 
    {
        auto _decompressed_stream = getDecompressedAlignmentStream();
        auto alignments_with_offsets = alignmentRange!withOffsets(_decompressed_stream);

        static struct Result(alias IteratePolicy, R, S) {
            this(R range, S stream, void delegate(lazy float p) progressBarFunc) {
                _range = range;
                _stream = stream;
                _progress_bar_func = progressBarFunc;
            }

            static if (__traits(identifier, IteratePolicy) == "withOffsets") {
                auto front() @property {
                    return _range.front;
                } 
            } else static if (__traits(identifier, IteratePolicy) == "withoutOffsets") {
                auto front() @property {
                    return _range.front.alignment;
                }
            } else static assert(0, __traits(identifier, IteratePolicy));

            bool empty() @property {
                return _range.empty;
            }

            void popFront() {
                _bytes_read += _range.front.alignment.size_in_bytes;
                _range.popFront();
                if (_progress_bar_func !is null) {
                    _progress_bar_func(min(1.0, 
                        cast(float)_bytes_read / (_stream.compressed_file_size * 
                                                  _stream.average_compression_ratio)));
                }
            }

            private R _range;
            private S _stream;
            private size_t _bytes_read;
            private void delegate(lazy float p) _progress_bar_func;
        }

        return Result!(IteratePolicy, 
                       typeof(alignments_with_offsets),
                       typeof(_decompressed_stream))(alignments_with_offsets, 
                                                     _decompressed_stream,
                                                     progressBarFunc);
    }

    /**
      Get an alignment at a given virtual offset.
     */
    Alignment getAlignmentAt(VirtualOffset offset) {
        enforce(_random_access_manager !is null);
        return _random_access_manager.getAlignmentAt(offset);
    }

    /**
      Get all alignments between two virtual offsets.

      First offset must point to the start of an alignment record,
      and be strictly less than the second one.

      For decompression, uses task pool specified at BamFile construction.
     */ 
    auto getAlignmentsBetween(VirtualOffset from, VirtualOffset to) {
        enforce(from < to, "First offset must be strictly less than second");
        enforce(_stream_is_seekable, "Stream is not seekable");
        
        return _random_access_manager.getAlignmentsBetween(from, to, _task_pool);
    }

    /**
      Get BAI chunks containing all reads overlapping specified region.
     */
    Chunk[] getChunks(int ref_id, int beg, int end) {
        enforce(_random_access_manager !is null);
        enforce(beg < end);

        return _random_access_manager.getChunks(ref_id, beg, end);
    }

    /**
      Returns reference sequence with id $(D ref_id).
     */
    ReferenceSequence reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return ReferenceSequence(_random_access_manager, 
                                 ref_id,
                                 _reference_sequences[ref_id]);
    }

    /**
      Returns reference sequence named $(D ref_name).
     */
    ReferenceSequence opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " does not exist");
        auto ref_id = _reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    /**
      Check if reference named $(D ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        return null != (ref_name in _reference_sequence_dict);
    }

    /**
      Set buffer size for I/O operations.
     */
    void setBufferSize(size_t buffer_size) {
        this.buffer_size = buffer_size;
    }

private:
    
    string _filename;                       // filename (if available)
    Stream _source_stream;                  // compressed
    IChunkInputStream _decompressed_stream; // decompressed
    Stream _bam;                            // decompressed + endian conversion
    bool _stream_is_seekable;

    // Virtual offset at which alignment records start.
    VirtualOffset _alignments_start_voffset;

    BaiFile _dont_access_me_directly_use_bai_file_for_that;
    enum BaiStatus {
        notInitialized,
        initialized,
        fileNotFound
    }
    BaiStatus _bai_status = BaiStatus.notInitialized;

    // provides access to index file
    @property ref BaiFile _bai_file() { // initialized lazily
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
        return _dont_access_me_directly_use_bai_file_for_that;
    }; 

    RandomAccessManager _rndaccssmgr; // unreadable for a purpose
    @property RandomAccessManager _random_access_manager() {
        if (_rndaccssmgr is null) {
            synchronized {
                auto bai = _bai_file; 
                // remember that it's lazily initialized,
                // so we need to do that to get the right BAI status

                if (_bai_status == BaiStatus.initialized) {
                    _rndaccssmgr = new RandomAccessManager(_filename, bai);
                } else {
                    _rndaccssmgr = new RandomAccessManager(_filename);
                }
            }
        }
        return _rndaccssmgr;
    }

    SamHeader _header;
    string _headertext; // for lazy SAM header parsing
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict; /// name -> index mapping

    TaskPool _task_pool;
    size_t buffer_size = 8192; // buffer size to be used for I/O

    Stream getNativeEndianSourceStream() {
        assert(_filename !is null);
        return new utils.stream.File(_filename);
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
    IChunkInputStream getDecompressedStream() {

        auto compressed_stream = getSeekableCompressedStream();

        auto bgzf_range = (compressed_stream is null) ? BgzfRange(_source_stream) :
                                                        BgzfRange(compressed_stream);
        version(serial) {
            auto chunk_range = map!decompressBgzfBlock(bgzf_range);
        } else {
            auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 24);
        }

        if (compressed_stream !is null) {
            return makeChunkInputStream(chunk_range, cast(size_t)compressed_stream.size);
        } else {
            return makeChunkInputStream(chunk_range);
        }
    }


    // get decompressed stream starting from the first alignment record
    IChunkInputStream getDecompressedAlignmentStream() {
        auto compressed_stream = getSeekableCompressedStream();

        if (compressed_stream !is null) {
            enforce(_alignments_start_voffset != 0UL);

            compressed_stream.seekCur(_alignments_start_voffset.coffset);
            auto bgzf_range = BgzfRange(compressed_stream);

            version(serial) {
                auto chunk_range = map!decompressBgzfBlock(bgzf_range);
            } else {
                auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 24);
            }
        
            auto sz = compressed_stream.size;
            auto stream = makeChunkInputStream(chunk_range, cast(size_t)sz);
            stream.readString(_alignments_start_voffset.uoffset);
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

        _headertext = to!string(_bam.readString(l_text));
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
