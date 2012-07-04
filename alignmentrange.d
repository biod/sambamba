module alignmentrange;

import virtualoffset;
import chunkinputstream;
import alignment;

import std.stream;
import std.algorithm;
import std.system;

struct RawAlignmentBlock {
    VirtualOffset start_virtual_offset;
    VirtualOffset end_virtual_offset;
    ubyte[] raw_alignment_data;
}

/// Whether to provide offsets or not
/// when iterating raw alignment blocks
private enum IteratePolicy {
    withOffsets,
    withoutOffsets
}

/**
  Range for iterating over unparsed alignments,
  together with their virtual offsets in the BAM file.
 */
private class UnparsedAlignmentRange(IteratePolicy policy) 
{ 

    /// Create new range from IChunkInputStream.
    this(ref IChunkInputStream stream) {
        _stream = stream;
        _endian_stream = new EndianStream(_stream, Endian.littleEndian);
        readNext();
    }

    bool empty() @property const {
        return _empty;
    }

static if (policy == IteratePolicy.withOffsets) {
    /**
        Returns: tuple of 1) virtual offset of alignment block beginning,
                          2) alignment block except the first 4 bytes (block_size)

        See SAM/BAM specification for details about the structure of alignment blocks.
     */
    RawAlignmentBlock front() @property {
        return RawAlignmentBlock(_start_voffset, 
                                 _stream.virtualTell(),
                                 _current_record);
    }
} else {
    /**
        Returns: alignment block except the first 4 bytes (block_size)
     */
    ubyte[] front() @property {
        return _current_record;
    }
}

    void popFront() {
        readNext();
    }

private:
    IChunkInputStream _stream;
    EndianStream _endian_stream;

    ubyte[] _current_record;
    bool _empty = false;

static if (policy == IteratePolicy.withOffsets) {
    VirtualOffset _start_voffset;
}

    /**
      Reads next alignment block from stream.
     */
    void readNext() {
        if (_stream.eof()) {
            _empty = true;
            return;
        }
        
static if (policy == IteratePolicy.withOffsets) {
        _start_voffset = _stream.virtualTell();
}

        int block_size = void;
        _endian_stream.read(block_size);
        _current_record = _stream.readSlice(block_size);
    }

}

/**
    Returns: range for iterating over alignment blocks
 */
auto unparsedAlignments(alias Policy)(ref IChunkInputStream stream) {
    return new UnparsedAlignmentRange!(Policy)(stream);
}

/// Returns: an alignment constructed out of given chunk of memory.
///
/// No parsing occurs, it is done lazily.
Alignment makeAlignment(ubyte[] chunk) {
    return Alignment(chunk);
}

/// Tuple of virtual offset of the alignment, and the alignment itself.
struct AlignmentBlock {
    VirtualOffset start_virtual_offset;
    VirtualOffset end_virtual_offset;
    Alignment alignment;
}

/// Returns: lazy range of Alignment structs constructed from a given stream.
auto alignmentRange(ref IChunkInputStream stream) {

    return map!makeAlignment(unparsedAlignments!(IteratePolicy.withoutOffsets)(stream));
}

/// Returns: lazy range of AlignmentBlock structs constructed from a given stream.
///
/// Provides start and end virtual offsets together with alignments.
auto alignmentRangeWithOffsets(ref IChunkInputStream stream) {
    static AlignmentBlock makeAlignmentBlock(RawAlignmentBlock block) {
        return AlignmentBlock(block.start_virtual_offset,
                              block.end_virtual_offset,
                              makeAlignment(block.raw_alignment_data));
    }

    return map!makeAlignmentBlock(unparsedAlignments!(IteratePolicy.withOffsets)(stream));
}
