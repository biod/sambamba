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
module alignmentrange;

import virtualoffset;
import chunkinputstream;
import alignment;

import std.stream;
import std.algorithm;
import std.system;
import utils.switchendianness;

/// Tuple of virtual offset of the alignment, and the alignment itself.
struct AlignmentBlock {
    VirtualOffset start_virtual_offset;
    VirtualOffset end_virtual_offset;
    Alignment alignment;
    alias alignment this;

    AlignmentBlock dup() @property const {
        return AlignmentBlock(start_virtual_offset, end_virtual_offset, alignment.dup);
    }
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

        // In fact, on BAM files containing a special EOF BGZF block
        // this condition will be always false!
        //
        // The reason is that we don't want to unpack next block just
        // in order to see if it's an EOF one or not.
        if (_stream.eof()) {
            _empty = true;
            return;
        }
     
        // In order to get the right virtual offset, we need to do it here.
        beforeNextAlignmentLoad();

        // Here's where _empty is really set!
        int block_size = void;
        ubyte* ptr = cast(ubyte*)(&block_size);
        auto _read = 0;
        while (_read < int.sizeof) {
            auto _actually_read = _endian_stream.readBlock(ptr, int.sizeof - _read);
            if (_actually_read == 0) {
                version(development) {
                    import std.stdio;
                    stderr.writeln("[info][alignment range] empty, read ", _read, " bytes, expected ", int.sizeof);
                }
                _empty = true;
                return;
            }
            _read += _actually_read;
            ptr += _actually_read;
        } 

        if (std.system.endian != Endian.littleEndian) {
            switchEndianness(&block_size, int.sizeof);
        }

        _current_record = Alignment(_stream.readSlice(block_size));
    }
}

/// Returns: lazy range of Alignment structs constructed from a given stream.
auto alignmentRange(alias IteratePolicy=withoutOffsets)(ref IChunkInputStream stream) {
    return new AlignmentRange!IteratePolicy(stream);
}
