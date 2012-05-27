/**
  Module for random access operations on BAM file.
 */
module randomaccessmanager;

import bgzfrange;
import virtualoffset;
import alignment;
import chunkinputstream;
import bai.read;
import bai.bin;
import bai.chunk;
import bai.utils.algo;

import std.stream;
import std.system;
import std.algorithm;
import std.array;
import std.range;
import std.traits;
import std.exception;

private {
    immutable int LINEAR_INDEX_WINDOW_SIZE_LOG = 14;
}

debug {
    import std.stdio;
}

// TODO: add caching
class RandomAccessManager {

    /// Constructs new manager for the stream.
    this(Stream stream) {
        _compressed_stream = new EndianStream(stream, Endian.littleEndian);
    }

    /// Constructs new manager with given index file.
    /// This allows to do random-access interval queries.
    this(Stream stream, ref BaiFile bai) {
        _compressed_stream = new EndianStream(stream, Endian.littleEndian);
        _bai = bai;
        _has_index_file = true;
    }

    /// Get single alignment at a given virtual offset.
    Alignment getAlignmentAt(VirtualOffset offset) {

        _compressed_stream.seekSet(cast(size_t)(offset.coffset));

        auto bgzf_range = new BgzfRange(_compressed_stream);
        auto decompressed_range = map!decompressBgzfBlock(bgzf_range);

        IChunkInputStream stream = makeChunkInputStream(decompressed_range);
        stream.readString(offset.uoffset); // TODO: optimize

        return alignmentRange(stream).front;
    }

    bool has_index_file() @property {
        return _has_index_file;
    }
    private bool _has_index_file = false; // overwritten in constructor if filename is provided

    static struct RefSeq {

        RandomAccessManager manager;
        uint ref_id;

        this(RandomAccessManager manager, uint ref_id) {
            this.manager = manager;
            this.ref_id = ref_id;
        }

        auto getAlignments(int start, int end) {
            return manager.getAlignments(ref_id, start, end);
        }
    }

    /// Returns an object representing reference sequence.
    auto reference_sequence(uint ref_id) {

        return RefSeq(this, ref_id);
    }

    /// Fetch alignments with given reference sequence id, overlapping [beg..end)
    auto getAlignments(size_t ref_id, int beg, int end) {

        enforce(has_index_file, "BAM index file (.bai) must be provided");
        enforce(_bai.indices.length > ref_id, "Invalid reference sequence index");

        auto bgzf_range = new BgzfRange(_compressed_stream);
        auto decompressed_range = map!decompressBgzfBlock(bgzf_range);

        IChunkInputStream stream = makeChunkInputStream(decompressed_range);

        bool canOverlap(Bin b) {
            return b.canOverlapWith(beg, end);
        }

        //auto bins = filter!((Bin b) { return b.canOverlapWith(beg, end); })
        //                   (_bai.indices[ref_id].bins);
        auto bins = filter!canOverlap(_bai.indices[ref_id].bins);

        debug {
            writeln("Ref. ID: ", ref_id);
            writeln(beg, " .. ", end);
            write("Bins: ");
            auto _bins = filter!((Bin b) { return b.canOverlapWith(beg, end); })
                           (_bai.indices[ref_id].bins);
            foreach (bin; _bins) {
                write(bin.id, " ");
            }
            writeln();
        }

        auto chunks = joiner(map!"a.chunks"(bins));

        debug {
            writeln("Chunks:");
            auto _chunks = chunks.save;
            foreach (chunk; joiner(map!"a.chunks"(bins))) {
                auto b = VirtualOffset(chunk.beg);
                auto e = VirtualOffset(chunk.end);
                writeln("\t", b.coffset, "/", b.uoffset, " .. ", 
                              e.coffset, "/", e.uoffset);
            }
            writeln();
        }

        auto _beg = max(0, beg - 1); // 1-based => subtract one
        auto _i = min(_beg >> LINEAR_INDEX_WINDOW_SIZE_LOG, 
                      _bai.indices[ref_id].ioffsets.length - 1);
        auto min_offset = (_i == -1) ? -1UL : _bai.indices[ref_id].ioffsets[_i];

        debug {
            auto v = VirtualOffset(min_offset);
            writeln("Min. offset: ", v.coffset, "/", v.uoffset);
        }

        bool minOffsetFilter(Chunk a) {
            return a.end > min_offset;
        }
        auto filtered_chunks = filter!minOffsetFilter(chunks);
//        auto filtered_chunks = filter!((Chunk a){ return a.end > min_offset; })(chunks); 

        debug {
            writeln("Filtered chunks: ");
            auto _filtered_chunks = filter!((Chunk a){ return a.end > min_offset; })(chunks); 
            foreach (chunk; _filtered_chunks) {
                auto b = VirtualOffset(chunk.beg);
                auto e = VirtualOffset(chunk.end);
                writeln("\t", b.coffset, "/", b.uoffset, " .. ", 
                              e.coffset, "/", e.uoffset);
            }
        }
        
        auto chunks_array = array(filtered_chunks);
        sort(chunks_array);

        debug {
            writeln("Sorted chunks: ");
            foreach (chunk; chunks_array) {
                auto b = VirtualOffset(chunk.beg);
                auto e = VirtualOffset(chunk.end);
                writeln("\t", b.coffset, "/", b.uoffset, " .. ", 
                              e.coffset, "/", e.uoffset);
            }
            stdout.flush();
        }

        auto disjoint_chunks = nonOverlappingChunks(chunks_array);

        debug {
            writeln("Optimized chunks: ");
            auto _optimized = nonOverlappingChunks(chunks_array);
            foreach (chunk; _optimized) {
                auto b = VirtualOffset(chunk.beg);
                auto e = VirtualOffset(chunk.end);
                writeln("\t", b.coffset, "/", b.uoffset, " .. ", 
                              e.coffset, "/", e.uoffset);
            }
            stdout.flush();
        }

        auto alignments = disjointChunkAlignmentRange(disjoint_chunks);
        while (!alignments.empty() && alignments.front.ref_id < ref_id) {
            alignments.popFront();
        }

        return filterAlignments(alignments, beg, end);
    }

private:

    Stream _compressed_stream;
    BaiFile _bai;

    static struct DisjointChunkAlignmentRange(Range) {

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
                if (_current_alignment_block.virtual_offset >= _current_chunk.end) {
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

        /// (VirtualOffset, Alignment) tuples
        ReturnType!alignmentRangeWithOffsets _alignment_blocks;

        AlignmentBlock _current_alignment_block;

        /// setup new alignment range, starting from a given offset
        void setupStream(VirtualOffset offset) {

            /// seek coffset in compressed stream
            _compressed_stream.seekSet(cast(size_t)(offset.coffset));

            /// setup BgzfRange and ChunkInputStream
            auto bgzf_range = new BgzfRange(_compressed_stream);
            // TODO: decompressing should use caching
            auto decompressed_range = map!decompressBgzfBlock(bgzf_range);
            IChunkInputStream stream = makeChunkInputStream(decompressed_range);

            /// seek uoffset in decompressed stream
            stream.readString(offset.uoffset); // TODO: optimize

            /// Setup new alignment range, with offsets.
            /// Offsets are needed to stop when we pass the end of 
            /// the current chunk, and skip to the next one.
            _alignment_blocks = alignmentRangeWithOffsets(stream);
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
                setupStream(VirtualOffset(_current_chunk.beg));
            }
        }
    }


    /// Range for iterating alignments contained in supplied intervals.
    /// 
    /// Modifies _compressed_stream during iteration.
    auto disjointChunkAlignmentRange(Range)(Range r) 
        if (is(ElementType!Range == Chunk))
    {
        return DisjointChunkAlignmentRange!Range(r, _compressed_stream);
    }

    static struct AlignmentFilter(R) {
        this(R r, int beg, int end) {
            _range = r;
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
                    _current_alignment.sequence_length < _beg) 
                {
                    /// ends before beginning of the region
                    ///  [-----------)
                    ///               [beg .......... end)
                    _range.popFront();
                } else {
                    return; /// _current_alignment overlaps the region
                }
            }
            _empty = true; /// seems like _range is empty...
        }
    }

    /// Get range of alignments sorted by leftmost coordinate,
    /// together with an interval [beg, end),
    /// and return another range of alignments which overlap the region.
    auto filterAlignments(R)(R r, int beg, int end) 
        if(is(ElementType!R == Alignment)) 
    {
        return AlignmentFilter!R(r, beg, end);
    }

}

