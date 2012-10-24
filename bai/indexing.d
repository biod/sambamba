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
module bai.indexing;

import bamfile;
import constants;

import bai.bin;
import bai.chunk;

import std.stream;
import std.array;
import std.algorithm;
import std.system;
import std.exception;

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

void defaultProgressBarFunc(lazy float dummy) {}

/// Writes BAM index to the $(D stream)
///
/// Accepts optional $(D progressBarFunc)
void createIndex(ref BamFile bam, ref Stream stream, void delegate(lazy float p) progressBarFunc=null) {

    auto endian_stream = new EndianStream(stream, Endian.littleEndian);

    auto refs = bam.reference_sequences;
    auto nrefs = refs.length;

    endian_stream.writeString(BAI_MAGIC);            // write BAI magic string
    endian_stream.write(cast(int)nrefs);             // and number of references

    void writeEmptyReference() {
        endian_stream.write(cast(int)0); // n_bins
        endian_stream.write(cast(int)0); // n_intv
    }

    // BAM file contains no alignments at all or all reads are unmapped
    if (bam.alignments!withOffsets.empty ||
        bam.alignments!withOffsets.front.alignment.ref_id < 0) {
        foreach (i; 0 .. nrefs) {
            writeEmptyReference();
        }
        return;
    }

    // OK, now let's deal with non-degenerate case

    auto alignment_blocks = bam.alignmentsWithProgress!withOffsets(progressBarFunc);

    auto prev_block = alignment_blocks.front;
    alignment_blocks.popFront();

    // this is the main character hereafter
    auto prev_read = prev_block.alignment;

    // array of linear offsets for the current reference entry
    ulong[BAI_MAX_BIN_ID - BAI_MAX_NONLEAF_BIN_ID + 1] linear_index;
    // (maximum index in linear_index where data was written) + 1
    size_t linear_index_write_length;

    // map: bin ID -> array of chunks
    Chunk[][uint] chunks;

    auto first_ref_id = prev_block.alignment.ref_id;
    auto current_chunk_beg = prev_block.start_virtual_offset;
    assert(first_ref_id >= 0);

    foreach (i; 0 .. first_ref_id) {
        writeEmptyReference();
    }

    void updateLinearIndex() {
        assert(prev_read.ref_id >= 0);

        size_t beg, end;

        if (prev_read.is_unmapped) {
            end = beg = toLinearIndexOffset(prev_read.position);
        } else {
            beg = toLinearIndexOffset(prev_read.position);
            end = toLinearIndexOffset(prev_read.position + prev_read.basesCovered() - 1);
        }

        debug {
            import std.stdio;
            if (end >= linear_index.length) {
                writeln("beg: ", beg);
                writeln("end: ", end);
                writeln("pos: ", prev_read.position);
                writeln("bases: ", prev_read.basesCovered());
            }
        }

        foreach (i; beg .. end + 1) {
            if (linear_index[i] == 0UL) {
                linear_index[i] = cast(ulong)prev_block.start_virtual_offset;
            }
        }

        if (end + 1 > linear_index_write_length) {
            linear_index_write_length = end + 1;
        }
    }

    void dumpCurrentLinearIndex() {
        endian_stream.write(cast(int)linear_index_write_length);

        //                                                                 
        // There might be untouched places in linear index                 
        // with virtual offset equal to zero.                              
        // However, it's not a good idea to leave those zeros,             
        // since we can start lookup from the last non-zero virtual offset 
        // encountered before the untouched window.                        
        //                                                                 
        ulong last_voffset = 0;

        foreach (voffset; linear_index[0 .. linear_index_write_length])
        {
            if (voffset == 0) {
                voffset = last_voffset;
            } else {
                last_voffset = voffset;
            }
            endian_stream.write(voffset);
        }
    }

    void dumpCurrentReference() {
        endian_stream.write(cast(int)chunks.length);

        foreach (bin_id, bin_chunks; chunks) {
            if (bin_chunks.length > 0) {
                endian_stream.write(bin_id);
                endian_stream.write(cast(int)bin_chunks.length);
                foreach (chunk; bin_chunks) {
                    endian_stream.write(cast(ulong)chunk.beg);
                    endian_stream.write(cast(ulong)chunk.end);
                }
            }
        }

        dumpCurrentLinearIndex();

        // reset data
        linear_index[] = 0;
        linear_index_write_length = 0;
        chunks = null;
        current_chunk_beg = prev_block.end_virtual_offset;
    }

    // adds chunk to the current bin (which is determined from prev_read)
    void updateChunks() {
        auto current_chunk_end = prev_block.end_virtual_offset;

        auto bin_id = prev_read.bin.id;

        if (bin_id !in chunks) {
            chunks[bin_id] = [];
        }
        auto cs = chunks[bin_id];

        bool canMergeWithPreviousChunk() {
            assert(cs.length > 0);
            auto last_chunk = cs[$ - 1];

            if (last_chunk.end.coffset == current_chunk_beg.coffset)
                return true;

            return false;
        }

        if (cs.length == 0 || !canMergeWithPreviousChunk()) {
            chunks[prev_read.bin.id] ~= Chunk(current_chunk_beg, current_chunk_end);
        } else {
            chunks[prev_read.bin.id][$ - 1].end = current_chunk_end;
        }

        current_chunk_beg = current_chunk_end;
    }

    foreach (block; alignment_blocks) {

        auto read = block.alignment;

        // new reference, so write data for previous one(s)
        if (read.ref_id != prev_read.ref_id) {
            updateLinearIndex();
            updateChunks();
            dumpCurrentReference();
            
            foreach (i; prev_read.ref_id + 1 .. read.ref_id)
                writeEmptyReference();
        }

        // this and all the following reads are unmapped
        if (read.ref_id < 0) {
            break;
        }

        // start position is unavailable, skip
        if (read.position < 0) {
            prev_block = block;
            prev_read = read;
            continue;
        }

        // check if the BAM file is indeed sorted
        if ((read.ref_id == prev_read.ref_id && 
             read.position < prev_read.position) ||
            (read.ref_id < prev_read.ref_id)) 
        {
            throw new Exception("BAM file is not properly sorted: " ~
                                "read '" ~ read.read_name ~ "'" ~
                                " must be before read '" ~ 
                                prev_read.read_name ~ 
                                "' (at virtual offset " ~ 
                                to!string(prev_block.start_virtual_offset)~
                                ")");
        }

        // ---------------------------------------------------------------------

        if (read.ref_id == prev_read.ref_id) {
            updateLinearIndex();

            if (read.bin.id != prev_read.bin.id) {
                updateChunks();
            }
        }

        // ---------------------------------------------------------------------

        prev_block = block;
        prev_read = read;
    }

    // after the loop, prev_read is the last read with ref_id >= 0
    assert(prev_read.ref_id >= 0);
    updateLinearIndex();
    updateChunks();
    dumpCurrentReference();

    // write the rest
    foreach (i; prev_read.ref_id + 1 .. nrefs) {
        writeEmptyReference();
    }
}
