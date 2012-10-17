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
module reference;

import randomaccessmanager;
import alignmentrange;
import virtualoffset;

import std.stream;
import std.exception;
import std.array;

/**
  Stores reference sequence name and length
 */
struct ReferenceSequenceInfo {
    string name;
    int length;

    /**
      Constructs the structure from input stream
     */
    this(ref Stream stream) {
        int l_name; // length of the reference name plus one
        stream.read(l_name);
        name = stream.readString(l_name)[0..$-1].idup; // strip '\0' at the end
        stream.read(length);
    }
}

/**
  Represents reference sequence.
 */
struct ReferenceSequence {
   
    /// Name of reference sequence as in BAM file
    string name() @property const {
        return _info.name;
    }

    /// Length in base pairs
    int length() @property const {
        return _info.length;
    }

    /// Reference ID
    int id() @property const {
        return _ref_id;
    }

    /// Get alignments overlapping [start, end)
    auto opSlice(int start, int end) {
        enforce(start < end, "start must be less than end");
        enforce(_manager !is null, "random access is not available");
        return _manager.getAlignments(_ref_id, start, end);
    }

    /// Get all alignments
    auto opSlice() {
        return opSlice(0, length);
    }

    private alias typeof(opSlice().front) Read;
    private Read _first_read() @property {
        return opSlice().front.dup;
    }

    /// First position on the reference overlapped by reads (0-based)
    int firstPosition() {
        return _first_read.position;
    }

    /// Virtual offset at which reads aligned to this reference start.
    VirtualOffset startVirtualOffset() {
        return _first_read.start_virtual_offset;
    }

    /// Last position on the reference overlapped by reads (0-based)
    int lastPosition() {
        // The key idea is
        //  1) use last offset from linear index
        //  2) loop through all remaining reads starting from there

        auto ioffsets = _manager.getBai().indices[_ref_id].ioffsets[];

        long index = ioffsets.length - 1;

        debug {
            int reads_processed = 0;
        }

        while (index >= 0) {
            auto offset = ioffsets[index];

            auto stream = _manager.createStreamStartingFrom(offset);
            auto reads = alignmentRange(stream);

            int last_position = int.min;

            foreach (read; reads) {

                 debug {
                     reads_processed += 1;
                 }

                 if (read.ref_id != _ref_id) {
                     break;
                 }
                
                 if (read.position == -1) {
                     continue;
                 }

                 auto end_pos = read.position + read.basesCovered();
                 if (end_pos > last_position)
                     last_position = end_pos;
            }

            if (last_position != int.min) {
                debug {
                    import std.stdio;
                    stderr.writeln("[debug] ReferenceSequence.lastPosition() processed ",
                                   reads_processed, " reads");
                }
                return last_position - 1;
            }

            --index;
        }

        return firstPosition();
    }

    this(RandomAccessManager manager, int ref_id, ReferenceSequenceInfo info) {
        _manager = manager;
        _ref_id = ref_id;
        _info = info;
    }

private:
    RandomAccessManager _manager;
    int _ref_id;
    ReferenceSequenceInfo _info;
}
