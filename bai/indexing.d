module bai.indexing;

import bamfile;
import constants;

import bai.bin;
import bai.chunk;

import std.stream;
import std.array;
import std.system;
import std.exception;

private size_t toLinearIndexOffset(int position) {
    return position == 0 ? 0 : (position - 1) >> BAI_LINEAR_INDEX_SHIFT;
}

/// Writes BAM index to the $(D stream)
void createIndex(ref BamFile bam, ref Stream stream) {

    auto endian_stream = new EndianStream(stream, Endian.littleEndian);

    auto refs = bam.reference_sequences;
    auto nrefs = refs.length;

    endian_stream.writeString(BAI_MAGIC);            // write BAI magic string
    endian_stream.write(cast(int)nrefs);       // and number of references

    void writeEmptyReference() {
        endian_stream.write(cast(int)0); // n_bins
        endian_stream.write(cast(int)0); // n_intv
    }

    // BAM file contains no alignments at all or all reads are unmapped
    if (bam.alignmentsWithVirtualOffsets.empty ||
        bam.alignmentsWithVirtualOffsets.front.alignment.ref_id < 0) {
        foreach (i; 0 .. nrefs) {
            writeEmptyReference();
        }
        return;
    }

    // OK, now let's deal with non-degenerate case

    auto alignment_blocks = bam.alignmentsWithVirtualOffsets;

    auto prev_block = alignment_blocks.front;
    alignment_blocks.popFront();

    // this is the main character hereafter
    auto prev_read = prev_block.alignment;

    // array of linear offsets for the current reference entry
    ulong[BAI_MAX_NONLEAF_BIN_ID + 1] linear_offsets;
    // (maximum index in linear_offsets where data was written) + 1
    size_t linear_offsets_write_length;

    // map: bin ID -> array of chunks
    Chunk[][uint] chunks;

    auto first_ref_id = prev_block.alignment.ref_id;
    auto current_chunk_beg = prev_block.start_virtual_offset;
    assert(first_ref_id >= 0);

    foreach (i; 0 .. first_ref_id) {
        writeEmptyReference();
    }

    void updateLinearOffsets() {
        assert(prev_read.ref_id >= 0);

        if (!prev_read.bin.is_leaf) {

            auto beg = toLinearIndexOffset(prev_read.position);
            auto end = toLinearIndexOffset(prev_read.position + prev_read.basesCovered());
            // TODO: think about moving basesCovered() calculation in another thread

            foreach (i; beg .. end + 1) {
                if (linear_offsets[i] == 0UL) {
                    linear_offsets[i] = cast(ulong)prev_block.start_virtual_offset;
                }
            }

            if (end + 1 > linear_offsets_write_length) {
                linear_offsets_write_length = end + 1;
            }
        }
    }

    void dumpCurrentLinearOffsets() {
        endian_stream.write(cast(int)linear_offsets_write_length);
        foreach (voffset; linear_offsets[0 .. linear_offsets_write_length])
        {
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

        dumpCurrentLinearOffsets();

        // reset data
        linear_offsets[] = 0;
        linear_offsets_write_length = 0;
        chunks = null;
        current_chunk_beg = prev_block.end_virtual_offset;
    }

    // adds chunk to the current bin (which is determined from prev_read)
    void updateChunks() {
        auto current_chunk_end = prev_block.end_virtual_offset;

        // add chunk
        chunks[prev_read.bin.id] ~= Chunk(current_chunk_beg, current_chunk_end);

        current_chunk_beg = current_chunk_end;
    }

    foreach (block; alignment_blocks) {

        auto read = block.alignment;

        // new reference, so write data for previous one(s)
        if (read.ref_id != prev_read.ref_id) {
            updateLinearOffsets();
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
                                to!string(prev_block.start_virtual_offset));
        }

        // ---------------------------------------------------------------------

        updateLinearOffsets();

        if (read.bin.id != prev_read.bin.id) {
            updateChunks();
        }

        // ---------------------------------------------------------------------

        prev_block = block;
        prev_read = read;
    }

    // after the loop, prev_read is the last read with ref_id >= 0
    assert(prev_read.ref_id >= 0);
    updateLinearOffsets();
    updateChunks();
    dumpCurrentReference();

    // write the rest
    foreach (i; prev_read.ref_id + 1 .. nrefs) {
        writeEmptyReference();
    }
}
