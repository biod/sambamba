/*
    This file is part of BioD.
    Copyright (C) 2012-2017    Artem Tarasov <lomereiter@gmail.com>

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
module bio.std.hts.bam.bai.indexing;

import bio.std.hts.bam.reader;
import bio.std.hts.bam.readrange;
import bio.std.hts.bam.constants;

import bio.std.hts.bam.bai.bin;
import bio.core.bgzf.chunk;

import contrib.undead.stream;
import std.array;
import std.algorithm;
import std.system;
import std.exception;
import core.stdc.string;

// Suppose we have an alignment which covers bases on a reference,
// starting from one position and ending at another position.
// In order to build linear index, we need to find to which windows
// the two positions correspond. 
//
//
// (K = 16384)
//
// [0, K)[K, 2K)[2K, 3K)...         <- windows
//    [.......)                     <- alignment
// 
private size_t toLinearIndexOffset(int position) {
    return position < 0 ? 0 : position / BAI_LINEAR_INDEX_WINDOW_SIZE;
}

///
struct IndexBuilder {
    private {
        // array of linear offsets for the current reference entry
        VirtualOffset[] _linear_index;
        // (maximum index in _linear_index where data was written) + 1
        size_t _linear_index_write_length;

        static struct PreviousRead {
            int ref_id = -1;
            int position;
            int end_position;
            int basesCovered() { return end_position - position; }
            Bin bin;
            bool is_unmapped;

            private char[256] _name_buf; // by spec., name length is <= 255 chars
            string name;

            VirtualOffset start_virtual_offset;
            VirtualOffset end_virtual_offset;
        }

        PreviousRead _prev_read;

        Stream _stream;
        int _n_refs;

        ulong _no_coord = 0;
        // metadata for each reference
        ulong _beg_vo = -1UL;
        ulong _end_vo = 0;
        ulong _unmapped = 0;
        ulong _mapped = 0;

        bool _first_read = true;

        // map: bin ID -> array of chunks
        Chunk[][uint] _chunks;

        VirtualOffset _current_chunk_beg; // start of current chunk

        // no metadata for empty references
        void writeEmptyReference() { 
            _stream.write(cast(int)0); // n_bins
            _stream.write(cast(int)0); // n_intv
        }

        void updateLastReadInfo(ref BamReadBlock read) {
            with (_prev_read) {
                ref_id = read.ref_id;
                position = read.position;
                end_position = position + read.basesCovered();
                bin = read.bin;
                is_unmapped = read.is_unmapped;
                _name_buf[0 .. read.name.length] = read.name[];
                name = cast(string)_name_buf[0 .. read.name.length];
                start_virtual_offset = read.start_virtual_offset;
                end_virtual_offset = read.end_virtual_offset;
            }
        }

        void updateMetadata(ref BamReadBlock read) {
            if (read.ref_id == -1) {
                ++_no_coord;
            } else {
                if (read.is_unmapped) {
                    ++_unmapped;
                } else {
                    ++_mapped;
                }

                if (_beg_vo == -1UL)
                    _beg_vo = cast(ulong)read.start_virtual_offset;
                _end_vo = cast(ulong)read.end_virtual_offset;
            }
        }

        void updateLinearIndex() {
            assert(_prev_read.ref_id >= 0);

            size_t beg, end;

            if (_prev_read.is_unmapped) {
                end = beg = toLinearIndexOffset(_prev_read.position);
            } else {
                beg = toLinearIndexOffset(_prev_read.position);
                end = toLinearIndexOffset(_prev_read.position + _prev_read.basesCovered() - 1);
            }

            debug {
                import std.stdio;
                if (end >= _linear_index.length) {
                    writeln("beg: ", beg);
                    writeln("end: ", end);
                    writeln("pos: ", _prev_read.position);
                    writeln("bases: ", _prev_read.basesCovered());
                }
            }

            foreach (i; beg .. end + 1)
                if (_linear_index[i] == 0UL)
                    _linear_index[i] = _prev_read.start_virtual_offset;

            if (end + 1 > _linear_index_write_length)
                _linear_index_write_length = end + 1;
        }

        void dumpCurrentLinearIndex() {
            _stream.write(cast(int)_linear_index_write_length);

            //                                                                 
            // There might be untouched places in linear index                 
            // with virtual offset equal to zero.                              
            // However, it's not a good idea to leave those zeros,             
            // since we can start lookup from the last non-zero virtual offset 
            // encountered before the untouched window.                        
            //                                                                 
            VirtualOffset last_voffset = 0;

            foreach (voffset; _linear_index[0 .. _linear_index_write_length]) {
                if (voffset == 0)
                    voffset = last_voffset;
                else
                    last_voffset = voffset;
                _stream.write(cast(ulong)voffset);
            }
        }

        void dumpCurrentReference() {
            // +1 because we output dummy bin, too
            _stream.write(cast(int)(_chunks.length + 1));

            foreach (bin_id, bin_chunks; _chunks) {
                if (bin_chunks.length > 0) {
                    _stream.write(cast(uint)bin_id);
                    _stream.write(cast(int)bin_chunks.length);
                    foreach (chunk; bin_chunks) {
                        _stream.write(cast(ulong)chunk.beg);
                        _stream.write(cast(ulong)chunk.end);
                    }
                }
            }
            _stream.write(cast(uint)37450);
            _stream.write(cast(int)2);
            _stream.write(cast(ulong)_beg_vo);
            _stream.write(cast(ulong)_end_vo);
            _stream.write(cast(ulong)_mapped);
            _stream.write(cast(ulong)_unmapped);

            dumpCurrentLinearIndex();

            // reset data
            memset(_linear_index.ptr, 0, _linear_index.length * ulong.sizeof);
            _linear_index_write_length = 0;
            _chunks = null;
            _current_chunk_beg = _prev_read.end_virtual_offset;

            _beg_vo = _end_vo = cast(ulong)_current_chunk_beg;
            _unmapped = 0;
            _mapped = 0;
        }

        // adds chunk to the current bin (which is determined from _prev_read)
        void updateChunks() {
            auto current_chunk_end = _prev_read.end_virtual_offset;

            auto bin_id = _prev_read.bin.id;

            if (bin_id !in _chunks)
                _chunks[bin_id] = [];
            auto cs = _chunks[bin_id];

            bool canMergeWithPreviousChunk() {
                assert(cs.length > 0);
                auto last_chunk = cs[$ - 1];

                if (last_chunk.end.coffset == _current_chunk_beg.coffset)
                    return true;

                return false;
            }

            if (cs.length == 0 || !canMergeWithPreviousChunk()) {
                auto new_chunk = Chunk(_current_chunk_beg, current_chunk_end);
                _chunks[_prev_read.bin.id] ~= new_chunk;
            } else {
                _chunks[_prev_read.bin.id][$ - 1].end = current_chunk_end;
            }

            _current_chunk_beg = current_chunk_end;
        }

        void checkThatBinIsCorrect(ref BamReadBlock read) {
            if (!check_bins)
                return;
            auto expected = reg2bin(read.position, 
                                    read.position + read.basesCovered());
            enforce(read.bin.id == expected,
                    "Bin in read with name '" ~ read.name ~ 
                    "' is set incorrectly (" ~ to!string(read.bin.id) ~ 
                    " instead of expected " ~ to!string(expected) ~ ")");
        }

        void checkThatInputIsSorted(ref BamReadBlock read) {
            if (_first_read) return;
            if (read.ref_id == -1) return; // unmapped
            if (_prev_read.ref_id < read.ref_id) return;

            enforce(read.ref_id == _prev_read.ref_id && 
                    read.position >= _prev_read.position,
                    "BAM file is not coordinate-sorted: " ~
                    "read '" ~ read.name ~ "' (" ~ read.ref_id.to!string ~ ":" ~ read.position.to!string ~ ")" ~
                    " must be after read '" ~ _prev_read.name ~ "' (" ~ _prev_read.ref_id.to!string ~ ":" ~ _prev_read.position.to!string ~ ")" ~
                    "' (at virtual offsets " ~ 
                    to!string(_prev_read.start_virtual_offset) ~ ", " ~ read.start_virtual_offset.to!string ~ ")");
        }
    }

    ///
    this(Stream output_stream, int number_of_references) {
        _stream = new EndianStream(output_stream, Endian.littleEndian);
        _n_refs = number_of_references;

        size_t size = BAI_MAX_BIN_ID - BAI_MAX_NONLEAF_BIN_ID + 1;
        _linear_index = new VirtualOffset[](size);

        _stream.writeString(BAI_MAGIC);            // write BAI magic string
        _stream.write(cast(int)_n_refs);           // and number of references
    }

    /// Check that bins are correct.
    bool check_bins;

    /// Add a read. The reads must be put in coordinate-sorted order.
    void put(BamReadBlock read) {
        checkThatInputIsSorted(read);
        scope(exit) updateMetadata(read);

        if (read.ref_id < 0)
            return;

        // start position is unavailable, skip
        if (read.position < 0)
            return;

        if (_first_read) {
            updateLastReadInfo(read);
            _first_read = false;
            _current_chunk_beg = read.start_virtual_offset;

            if (read.ref_id > 0)
                foreach (i; 0 .. read.ref_id)
                    writeEmptyReference();

            return;
        }

        checkThatBinIsCorrect(read);

        // new reference, so write data for previous one(s)
        if (read.ref_id > _prev_read.ref_id) {
            updateLinearIndex();
            updateChunks();
            dumpCurrentReference();
            
            foreach (i; _prev_read.ref_id + 1 .. read.ref_id)
                writeEmptyReference();
        }

        if (read.ref_id == _prev_read.ref_id) {
            updateLinearIndex();

            if (read.bin.id != _prev_read.bin.id)
                updateChunks();
        }

        updateLastReadInfo(read);
    }

    /// Closes the stream
    void finish() {
        if (!_first_read) { // at least one was processed
            assert(_prev_read.ref_id >= 0);
            updateLinearIndex();
            updateChunks();
            dumpCurrentReference();
        }

        // _prev_read.ref_id == -1 if all are unmapped
        foreach (i; _prev_read.ref_id + 1 .. _n_refs)
            writeEmptyReference();

        _stream.write(cast(ulong)_no_coord);
        _stream.close();
    }
}

/// Writes BAM index to the $(D stream)
///
/// Accepts optional $(D progressBarFunc)
void createIndex(BamReader bam, Stream stream, bool check_bins=false,
                 void delegate(lazy float p) progressBarFunc=null)
{
    auto n_refs = cast(int)bam.reference_sequences.length;
    auto index_builder = IndexBuilder(stream, n_refs);
    index_builder.check_bins = check_bins;
    auto reads = bam.readsWithProgress!withOffsets(progressBarFunc);
    foreach (read; reads)
        index_builder.put(read);
    index_builder.finish();
}
