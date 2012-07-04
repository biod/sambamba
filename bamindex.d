import bamfile;
import constants;

import bai.bin;
import bai.chunk;

import std.stream;
import std.array;
import std.system;
import std.exception;

/// Writes BAM index to $(D stream)
///
/// Sets endianness of $(D stream) to Endian.littleEndian
void createIndex(ref BamFile bam, ref EndianStream stream) {

    stream.endian = Endian.littleEndian;

    auto refs = bam.reference_sequences;
    auto nrefs = refs.length;

    stream.write(BAI_MAGIC);            // write BAI magic string
    stream.write(cast(int)nrefs);       // and number of references

    void writeEmptyReference() {
        stream.write(cast(int)0); // n_bins
        stream.write(cast(int)0); // n_intv
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
    auto prev_read = prev_block.read;

    // array of linear offsets for the current reference entry
    uint[BAI_MAX_NONLEAF_BIN_ID + 1] linear_offsets;
    // (maximum index in linear_offsets where data was written) + 1
    size_t linear_offsets_write_length;

    // array of bins for the current reference entry
    Chunk[][uint] chunks;

    void updateLinearOffsets() {
        assert(prev_read.ref_id >= 0);

        if (!prev_read.bin.is_leaf) {

            auto beg = prev_read.position >> BAI_LINEAR_INDEX_SHIFT;
            auto end = (prev_read.position + prev_read.basesCovered() - 1) 
                            >> BAI_LINEAR_INDEX_SHIFT;

            // TODO: think about moving basesCovered() calculation in another thread

            for (i; beg .. end) {
                if (linear_offsets[i] == 0) {
                    linear_offsets[i] = prev_read.start_virtual_offset;
                }
            }
        }
    }

    void dumpBin(ref Bin bin) {
        stream.write(bin.id);
        stream.write(cast(int)bin.chunks.length);
        foreach (chunk; bin.chunks) {
            stream.write(chunk.beg);
            stream.write(chunk.end);
        }
    }

    void dumpCurrentReference() {
        stream.write(cast(int)bins.length);
        foreach (bin; bins) {
            if (bin.chunks.length > 0) {
                dumpBin(bin);
            }
        }

        dumpCurrentLinearOffsets();

        // reset data
        linear_offsets[] = 0;
        chunks = null;
    }

    auto first_ref_id = prev_block.alignment.ref_id;
    auto current_chunk_beg = prev_block.start_virtual_offset;
    assert(first_ref_id >= 0);

    foreach (i; 0 .. first_ref_id) {
        writeEmptyReference();
    }

    // adds chunk to the current bin (which is determined from prev_read)
    void updateChunks() {
        auto current_chunk_end = prev_read.end_virtual_offset;

        // add chunk
        chunks[prev_read.bin.id] ~= Chunk(current_chunk_beg, current_chunk_end);

        current_chunk_beg = current_chunk_end;
    }

    foreach (block; alignment_blocks) {

        auto read = block.alignment;

        // new reference, so write data for previous one(s)
        if (read.ref_id != prev_read.ref_id) {
            dumpCurrentReference();
            
            foreach (i; prev_read.ref_id + 1 .. read.ref_id)
                writeEmptyReference();
        }

        // this and all the following reads are unmapped
        if (read.ref_id < 0) {
            break;
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
