/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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
/**
  Module for random access operations on BAM file.
 */
module bio.std.hts.bam.randomaccessmanager;

import bio.std.hts.bam.constants;
import bio.std.hts.bam.reader;
import bio.std.hts.bam.read;
import bio.std.hts.bam.readrange;
import bio.std.hts.bam.baifile;
import bio.std.hts.bam.region;
import bio.core.utils.algo;

import bio.core.bgzf.block;
import bio.core.bgzf.virtualoffset;
import bio.core.bgzf.inputstream;
import bio.core.bgzf.constants;
import bio.core.bgzf.chunk;
import bio.core.utils.range;
import bio.core.utils.stream;

import std.system;
import std.algorithm;
import std.array;
import std.range;
import std.traits;
import std.exception;
import std.container;
import std.parallelism;
static import std.file;

debug {
    import std.stdio;
}

private {
    auto nonOverlappingChunks(R)(R chunks) {
        static ref auto chunkB(ref Chunk chunk) { return chunk.beg; }
        static ref auto chunkE(ref Chunk chunk) { return chunk.end; }
        return nonOverlapping!(chunkB, chunkE)(chunks);
    }
}

/// Class which random access tasks are delegated to.
class RandomAccessManager {

  // void setCache(BgzfBlockCache cache) {
      // _cache = cache;
  //}

    void setTaskPool(TaskPool task_pool) {
        _task_pool = task_pool;
    }

    void setBufferSize(size_t buffer_size) {
        _buffer_size = buffer_size;
    }

    /// Constructs new manager for BAM file
    this(string filename) {
        _filename = filename;
    }

    /// ditto
    this(BamReader reader) {
        _reader = reader;
        _filename = reader.filename;
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

    /// ditto
    this(BamReader reader, ref BaiFile bai) {
        _reader = reader;
        _filename = reader.filename;
        _bai = bai;
        _found_index_file = true;
    }

    /// If file ends with EOF block, returns virtual offset of the start of EOF block.
    /// Otherwise, returns virtual offset of the physical end of file.
    VirtualOffset eofVirtualOffset() const {
        ulong file_offset = std.file.getSize(_filename);
        if (hasEofBlock()) {
            return VirtualOffset(file_offset - BAM_EOF.length, 0);
        } else {
            return VirtualOffset(file_offset, 0);
        }
    }

    /// Returns true if the file ends with EOF block, and false otherwise.
    bool hasEofBlock() const {
        auto _stream = new bio.core.utils.stream.File(_filename);
        if (_stream.size < BAM_EOF.length) {
            return false;
        }

        ubyte[BAM_EOF.length] buf;
        _stream.seekEnd(-cast(int)BAM_EOF.length);

        _stream.readExact(&buf, BAM_EOF.length);
        if (buf != BAM_EOF) {
            return false;
        }

        return true;
    }

    /// Get new BgzfInputStream starting from specified virtual offset.
    BgzfInputStream createStreamStartingFrom(VirtualOffset offset)
    {
        auto _stream = new bio.core.utils.stream.File(_filename);
        auto _compressed_stream = new EndianStream(_stream, Endian.littleEndian);
        _compressed_stream.seekSet(cast(size_t)(offset.coffset));
        auto supplier = new StreamSupplier(_compressed_stream, offset.uoffset);
        // auto bgzf_stream = new BgzfInputStream(supplier, _task_pool, _cache);
        auto bgzf_stream = new BgzfInputStream(supplier, _task_pool);
        return bgzf_stream;
    }

    /// Get single read at a given virtual offset.
    /// Every time new stream is used.
    BamRead getReadAt(VirtualOffset offset) {
        auto stream = createStreamStartingFrom(offset);

        bool old_mode = _reader._seqprocmode;
        _reader._seqprocmode = true;
        auto read = bamReadRange(stream, _reader).front.dup;
        _reader._seqprocmode = old_mode;
        return read;
    }

    /// Get BGZF block at a given offset.
    BgzfBlock getBgzfBlockAt(ulong offset) {
        auto fstream = new bio.core.utils.stream.File(_filename);
        auto stream = new EndianStream(fstream, Endian.littleEndian);
        stream.seekSet(offset);
        BgzfBlock block = void;
        ubyte[BGZF_MAX_BLOCK_SIZE] buf = void;
        fillBgzfBufferFromStream(stream, true, &block, buf.ptr);
        block._buffer = block._buffer.dup;
        return block;
    }

    /// Get reads between two virtual offsets. First virtual offset must point
    /// to a start of an alignment record.
    auto getReadsBetween(VirtualOffset from, VirtualOffset to) {
        auto stream = createStreamStartingFrom(from);

        static bool offsetTooBig(BamReadBlock record, VirtualOffset vo) {
            return record.end_virtual_offset > vo;
        }

        return until!offsetTooBig(bamReadRange!withOffsets(stream, _reader), to);
    }

    bool found_index_file() @property const {
        return _found_index_file;
    }
    private bool _found_index_file = false; // overwritten in constructor if filename is provided

    /// BAI file
    ref const(BaiFile) getBai() const {
        enforce(found_index_file, "BAM index file (.bai) must be provided");
        return _bai;
    }

    private void checkIndexExistence() {
        enforce(found_index_file, "BAM index file (.bai) must be provided");
    }

    private void checkRefId(uint ref_id) {
        enforce(ref_id < _bai.indices.length, "Invalid reference sequence index");
    }

    private void appendChunks(ref Chunk[] chunks, Bin bin, VirtualOffset min_offset) {
        foreach (chunk; bin.chunks) {
            if (chunk.end > min_offset) {
                chunks ~= chunk;

                // optimization
                if (chunks[$-1].beg < min_offset)
                    chunks[$-1].beg = min_offset;
            }
        }
    }

    /// Get BAI chunks containing all alignment records overlapping the region
    Chunk[] getChunks(BamRegion region) {
        auto ref_id = region.ref_id;
        auto beg = region.start;
        auto end = region.end;
        checkIndexExistence();
        checkRefId(ref_id);

        // Select all bins that overlap with [beg, end).
        // Then from such bins select all chunks that end to the right of min_offset.
        // Sort these chunks by leftmost coordinate and remove all overlaps.

        auto min_offset = _bai.indices[ref_id].getMinimumOffset(beg);

        Chunk[] bai_chunks;
        foreach (b; _bai.indices[ref_id].bins) {
            if (!b.canOverlapWith(beg, end))
                continue;
            appendChunks(bai_chunks, b, min_offset);
        }

        sort(bai_chunks);
        return bai_chunks.nonOverlappingChunks().array();
    }

    // regions should be from the same reference sequence
    private Chunk[] getGroupChunks(BamRegion[] regions) {
        auto bitset = Array!bool();
        enforce(regions.length > 0);
        bitset.length = BAI_MAX_BIN_ID;
        bitset[0] = true;
        foreach (region; regions) {
            auto beg = region.start;
            auto end = region.end;
            int i = 0, k;
            enforce(beg < end);
            --end;
            for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) bitset[k] = true;
            for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) bitset[k] = true;
            for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) bitset[k] = true;
            for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) bitset[k] = true;
            for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) bitset[k] = true;
        }

        auto ref_id = regions.front.ref_id;
        checkIndexExistence();
        checkRefId(ref_id);

        // Select all bins that overlap with [beg, end).
        // Then from such bins select all chunks that end to the right of min_offset.
        // Sort these chunks by leftmost coordinate and remove all overlaps.

        auto min_offset = _bai.indices[ref_id].getMinimumOffset(regions.front.start);

        version(extraVerbose) {
            import std.stdio;
            stderr.writeln("min offset = ", min_offset);
        }

        Chunk[] bai_chunks;
        auto bins = _bai.indices[ref_id].bins;
        foreach (bin; bins)
            if (bitset[bin.id])
                appendChunks(bai_chunks, bin, min_offset);
        sort(bai_chunks);

        version(extraVerbose) {
            stderr.writeln("[chunks before normalization]");
            foreach(chunk; bai_chunks)
                stderr.writeln(chunk.beg, " - ", chunk.end);
        }

        return bai_chunks.nonOverlappingChunks().array();
    }

    private auto filteredReads(alias IteratePolicy)(BamRegion[] regions) {
        auto chunks = getGroupChunks(regions);
        version(extraVerbose) {
            import std.stdio;
            stderr.writeln("[chunks]");
            foreach (chunk; chunks)
                stderr.writeln(chunk.beg, " - ", chunk.end);
        }
        auto reads = readsFromChunks!IteratePolicy(chunks);
        return filterBamReads(reads, regions);
    }

    /// Fetch alignments with given reference sequence id, overlapping [beg..end)
    auto getReads(alias IteratePolicy=withOffsets)(BamRegion region)
    {
        auto chunks = getChunks(region);
        auto reads = readsFromChunks!IteratePolicy(chunks);
        return filterBamReads(reads, [region]);
    }

    auto getReads(alias IteratePolicy=withOffsets)(BamRegion[] regions) {
        auto sorted_regions = regions.sort();
        BamRegion[][] regions_by_ref;
        // TODO: replace with groupBy once it's included into Phobos
        uint last_ref_id = uint.max;
        foreach (region; sorted_regions) {
            if (region.ref_id == last_ref_id) {
                regions_by_ref.back ~= region;
            } else {
                regions_by_ref ~= [region];
                last_ref_id = region.ref_id;
            }
        }

        static ref auto regB(ref BamRegion region) { return region.start; }
        static ref auto regE(ref BamRegion region) { return region.end; }
        foreach (ref group; regions_by_ref)
            group = nonOverlapping!(regB, regE)(group).array();

        return regions_by_ref.zip(repeat(this))
                             .map!(gt => gt[1].filteredReads!IteratePolicy(gt[0]))()
                             .joiner();
    }

private:
    auto readsFromChunks(alias IteratePolicy, R)(R chunks) {
        auto fstream = new bio.core.utils.stream.File(_filename);
        auto compressed_stream = new EndianStream(fstream, Endian.littleEndian);
        auto supplier = new StreamChunksSupplier(compressed_stream, chunks);
        // auto stream = new BgzfInputStream(supplier, _task_pool, _cache);
        auto stream = new BgzfInputStream(supplier, _task_pool);
        return bamReadRange!IteratePolicy(stream, _reader);
    }

    string _filename;
    BaiFile _bai;
    BamReader _reader;
    TaskPool _task_pool;
    size_t _buffer_size;

  // BgzfBlockCache _cache;

    TaskPool task_pool() @property {
        if (_task_pool is null)
            _task_pool = taskPool;
        return _task_pool;
    }

public:

    static struct BamReadFilter(R) {
        this(R r, BamRegion[] regions) {
            _range = r;
            _regions = regions;
            enforce(regions.length > 0);
            _region = _regions.front;
            _ref_id = _region.ref_id; // assumed to be constant
            findNext();
        }

        bool empty() @property {
            return _empty;
        }

        ElementType!R front() @property {
            return _current_read;
        }

        void popFront() {
            _range.popFront();
            findNext();
        }

    private:
        R _range;
        uint _ref_id;
        BamRegion _region;
        BamRegion[] _regions; // non-overlapping and sorted
        bool _empty;
        ElementType!R _current_read;

        void findNext() {
            if (_regions.empty || _range.empty) {
                _empty = true;
                return;
            }

            while (!_range.empty) {
                _current_read = _range.front;

                // BamReads are sorted first by ref. ID.
                auto current_ref_id = cast(uint)_current_read.ref_id;
                // ref_id can't be -1 unless the index is fucked up
                if (current_ref_id > _ref_id) {
                    // no more records for this _ref_id
                    _empty = true;
                    return;
                } else if (current_ref_id < _ref_id) {
                    // skip reads referring to sequences
                    // with ID less than ours
                    _range.popFront();
                    continue;
                }

                if (_current_read.position >= _region.end) {
                    // As reads are sorted by leftmost coordinate,
                    // all remaining alignments in _range
                    // will not overlap the current interval as well.
                    //
                    //                  [-----)
                    //                  . [-----------)
                    //                  .  [---)
                    //                  .    [-------)
                    //                  .         [-)
                    //    [beg .....  end)
                    _regions.popFront();
                    // TODO: potentially binary search may be faster,
                    // but it needs to be checked
                    if (_regions.empty) {
                        _empty = true;
                        return;
                    } else {
                        _region = _regions.front;
                        continue;
                    }
                }

                if (_current_read.position > _region.start) {
                    return; // definitely overlaps
                }

                if (_current_read.position +
                    _current_read.basesCovered() <= _region.start)
                    {
                        /// ends before beginning of the region
                        ///  [-----------)
                        ///               [beg .......... end)
                        _range.popFront();
                        /// Zero-length reads are also considered non-overlapping,
                        /// so for consistency the inequality 12 lines above is strict.
                    } else {
                    return; /// _current_read overlaps the region
                }
            }
            _empty = true;
        }
    }

    // Get range of alignments sorted by leftmost coordinate,
    // together with an interval [beg, end),
    // and return another range of alignments which overlap the region.
    static auto filterBamReads(R)(R r, BamRegion[] regions)
    {
        return BamReadFilter!R(r, regions);
    }
}
