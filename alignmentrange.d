module alignmentrange;

import virtualoffset;
import chunkinputstream;
import alignment;

import std.stream;
import std.algorithm;
import std.system;

/// Tuple of virtual offset of the alignment, and the alignment itself.
struct AlignmentBlock {
    VirtualOffset start_virtual_offset;
    VirtualOffset end_virtual_offset;
    Alignment alignment;
}

/// Policies for alignmentRange
mixin template withOffsets() {
    /**
        Returns: virtual offsets of beginning and end of the current alignment
                 plus the current alignment itself.
     */
    AlignmentBlock front() @property {
        return AlignmentBlock(_start_voffset, 
                              _stream.virtualTell(),
                              _current_record);
    }

    private VirtualOffset _start_voffset;

    private void beforeNextAlignmentLoad() {
        _start_voffset = _stream.virtualTell();
    }
}

/// ditto
mixin template withoutOffsets() {
    /**
        Returns: current alignment
     */
    Alignment front() @property {
        return _current_record;
    }

    private void beforeNextAlignmentLoad() {}
}

class AlignmentRange(alias IteratePolicy) 
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

    mixin IteratePolicy;
    
    void popFront() {
        readNext();
    }

private:
    IChunkInputStream _stream;
    EndianStream _endian_stream;

    Alignment _current_record;
    bool _empty = false;

    /**
      Reads next alignment block from stream.
     */
    void readNext() {
        if (_stream.eof()) {
            _empty = true;
            return;
        }
       
        beforeNextAlignmentLoad();
        int block_size = void;
        _endian_stream.read(block_size);
        _current_record = Alignment(_stream.readSlice(block_size));
    }
}

/// Returns: lazy range of Alignment structs constructed from a given stream.
auto alignmentRange(alias IteratePolicy=withoutOffsets)(ref IChunkInputStream stream) {
    return new AlignmentRange!IteratePolicy(stream);
}
