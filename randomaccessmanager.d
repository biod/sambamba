/**
  Module for random access operations on BAM file.
 */
module randomaccessmanager;

import bgzfrange;
import virtualoffset;
import alignment;
import alignmentrange;
import chunkinputstream;
import bai.read;
import bai.bin;
import bai.chunk;
import bai.utils.algo;
import utils.memoize;
import utils.range;

import std.stream;
import std.system;
import std.algorithm;
import std.array;
import std.range;
import std.traits;
import std.exception;

import std.parallelism;

// made public to be accessible from utils.memoize module
auto decompressTask(BgzfBlock block) {
    auto t = task!decompressBgzfBlock(block);
    taskPool.put(t);
    return t;
}

private alias memoize!(decompressTask, 512, FifoCache, BgzfBlock) memDecompressTask;

// made public only to be accessible from std.algorithm
auto decompressSerial(BgzfBlock block) {
    return decompress(block).yieldForce();
}

// ditto
auto decompress(BgzfBlock block) { 
    return memDecompressTask(block);
}

debug {
    import std.stdio;
}

/// Class which random access tasks are delegated to.
class RandomAccessManager {

    /// Constructs new manager for BAM file
    this(string filename) {
        _filename = filename;
    }

    /// Constructs new manager with given index file.
    /// This allows to do random-access interval queries.
    ///
    /// Params:
    ///     filename =  location of BAM file
    ///     bai  =  index file
    this(string filename, ref BaiFile bai) {

        _filename = filename;
        _bai = bai;
        _found_index_file = true;
    }

    /// Get single alignment at a given virtual offset.
    /// Every time new stream is used.
    Alignment getAlignmentAt(VirtualOffset offset) {

        auto _stream = new BufferedFile(_filename);
        auto _compressed_stream = new EndianStream(_stream, Endian.littleEndian);
        _compressed_stream.seekSet(cast(size_t)(offset.coffset));

        auto bgzf_range = BgzfRange(_compressed_stream);
        auto decompressed_range = map!decompressSerial(bgzf_range);

        IChunkInputStream stream = makeChunkInputStream(decompressed_range);
        stream.readString(offset.uoffset); // TODO: optimize

        return alignmentRange(stream).front.dup;
    }

    bool found_index_file() @property {
        return _found_index_file;
    }
    private bool _found_index_file = false; // overwritten in constructor if filename is provided

    /// Fetch alignments with given reference sequence id, overlapping [beg..end)
    auto getAlignments(int ref_id, int beg, int end) {

        enforce(found_index_file, "BAM index file (.bai) must be provided");
        enforce(ref_id >= 0 && ref_id < _bai.indices.length, "Invalid reference sequence index");

        auto min_offset = _bai.indices[ref_id].getMinimumOffset(beg);

        auto _stream = new BufferedFile(_filename);
        Stream _compressed_stream = new EndianStream(_stream, Endian.littleEndian);

version(DigitalMars) {
        return _bai.indices[ref_id].bins
                                   .filter!((Bin b) { return b.canOverlapWith(beg, end); })
                                   .map!((Bin b) { return b.chunks; })
                                   .joiner()
                                   .filter!((Chunk c) { return c.end > min_offset; })
                                   .array()
                                   .sort()
                                   .nonOverlappingChunks()
                                   .disjointChunkAlignmentRange(_compressed_stream)
                                   .filterAlignments(ref_id, beg, end);

} else {
        /// use less functional approach
        Chunk[] chunks;
        foreach (b; _bai.indices[ref_id].bins) {
            if (!b.canOverlapWith(beg, end)) {
                continue;
            }

            foreach (chunk; b.chunks) {
                if (chunk.end > min_offset) {
                    chunks ~= chunk;
                }
            }
        }

        sort(chunks);

        auto disjoint_chunks = nonOverlappingChunks(chunks);

        auto alignments = disjointChunkAlignmentRange(disjoint_chunks, _compressed_stream);

        return filterAlignments(alignments, ref_id, beg, end);
} // endif

    }

private:
    
    string _filename;
    BaiFile _bai;
   
}

private {

    struct DisjointChunkAlignmentRange(Range) {

        this(Range r, ref Stream compressed_stream) {
            _chunks = r;
            _compressed_stream = compressed_stream;
            setupNextStream();
        }
        
        bool empty() @property {
            return _empty;
        }

        void popFront() {

            _alignment_blocks.popFront();

            if (!_alignment_blocks.empty) {
                _current_alignment_block = _alignment_blocks.front;

                /// Are we finished with the current chunk?
                if (_current_alignment_block.start_virtual_offset >= _current_chunk.end) {
                   _chunks.popFront();
                   setupNextStream();
                }

            } else {
                _empty = true; /// end of file
            }
        }

        Alignment front() @property {
            return _current_alignment_block.alignment;
        }

    private:
        Range _chunks;
        Chunk _current_chunk;
        Stream _compressed_stream;
        bool _empty = false;

        // start and end virtual offsets + alignments
        ReturnType!(alignmentRange!withOffsets) _alignment_blocks;

        AlignmentBlock _current_alignment_block;

        // setup new alignment range, starting from a given offset
        void setupStream(VirtualOffset offset) {

            // seek coffset in compressed stream
            _compressed_stream.seekSet(cast(size_t)(offset.coffset));

            // setup BgzfRange and ChunkInputStream
            auto bgzf_range = BgzfRange(_compressed_stream);

            version(serial) {
                auto decompressed_range = map!decompressBgzfBlock(bgzf_range);
            } else {
                /// up to 2 tasks are being executed at every moment
                auto prefetched_range = prefetch(map!decompress(bgzf_range), 2);
                auto decompressed_range = map!"a.yieldForce()"(prefetched_range);
            }
            IChunkInputStream stream = makeChunkInputStream(decompressed_range);

            /// seek uoffset in decompressed stream
            stream.readString(offset.uoffset); // TODO: optimize

            /// Setup new alignment range, with offsets.
            /// Offsets are needed to stop when we pass the end of 
            /// the current chunk, and skip to the next one.
            _alignment_blocks = alignmentRange!withOffsets(stream);
            if (!_alignment_blocks.empty) {
                _current_alignment_block = _alignment_blocks.front;

            } else {
                _empty = true;
            }
        }
      
        /// Tries to setup a stream corresponding to the front of chunk range.
        /// Side effects:
        ///         1) if _chunks is empty, sets _empty flag
        ///         2) otherwise, updates _current_chunk
        void setupNextStream() {
            if (_chunks.empty) {
                _empty = true;
            } else {
                _current_chunk = _chunks.front;
                setupStream(_current_chunk.beg);
            }
        }
    }


    // Range for iterating alignments contained in supplied intervals.
    // 
    // Modifies stream during iteration.
    auto disjointChunkAlignmentRange(Range)(Range r, ref Stream stream) 
        if (is(ElementType!Range == Chunk))
    {
        return DisjointChunkAlignmentRange!Range(r, stream);
    }
 
    struct AlignmentFilter(R) {
        this(R r, int ref_id, int beg, int end) {
            _range = r;
            _ref_id = ref_id;
            _beg = beg;
            _end = end;
            findNext();
        }

        bool empty() @property {
            return _empty;
        }

        Alignment front() @property {
            return _current_alignment;
        }
        
        void popFront() {
            _range.popFront();
            findNext();
        }

    private: 
        R _range;
        int _ref_id;
        int _beg;
        int _end;
        bool _empty;
        Alignment _current_alignment;

        void findNext() {
            if (_range.empty) {
                _empty = true;
                return;
            }
            while (!_range.empty) {
                _current_alignment = _range.front;

                /// Alignments are sorted first by ref. ID.
                auto current_ref_id = _current_alignment.ref_id;
                if (current_ref_id > _ref_id) {
                    /// no more records for this _ref_id
                    _empty = true;
                    return;
                } else if (current_ref_id < _ref_id) {
                    /// skip alignments referring to sequences
                    /// with ID less than ours
                    _range.popFront();
                    continue;
                }

                if (_current_alignment.position >= _end) {
                    _empty = true;
                    /// As alignments are sorted by leftmost
                    /// coordinate, all remaining alignments
                    /// in _range will not overlap the interval
                    /// as well.
                    /// 
                    ///                  [-----)
                    ///                  . [-----------)
                    ///                  .  [---)
                    ///                  .    [-------)
                    ///                  .         [-)
                    ///    [beg .....  end)
                    return;
                }

                if (_current_alignment.position +
                    _current_alignment.basesCovered() <= _beg) 
                {
                    /// ends before beginning of the region
                    ///  [-----------)
                    ///               [beg .......... end)
                    _range.popFront();
                } else {
                    return; /// _current_alignment overlaps the region
                }
            }
            _empty = true; 
        }
    }

    // Get range of alignments sorted by leftmost coordinate,
    // together with an interval [beg, end),
    // and return another range of alignments which overlap the region.
    auto filterAlignments(R)(R r, int ref_id, int beg, int end) 
        if(is(ElementType!R == Alignment)) 
    {
        return AlignmentFilter!R(r, ref_id, beg, end);
    }

}
