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
module bai.utils.algo;

import bai.chunk;

import std.range;
import std.algorithm;

struct NonOverlappingChunks(R) {

    this(R r) {
        _range = r;
        popFront();
    }

    bool empty() @property {
        return _empty;
    }

    auto front() @property {
        return _front;
    }

    void popFront() {
        if (!_range.empty()) {
            _front = _range.front;
            tryToJoinWithNextChunks();
        } else {
            _empty = true;
        }
    }

private:

    R _range;

    void tryToJoinWithNextChunks() {
        _range.popFront();
        while (!_range.empty()) {
            /// Get next element
            auto next = _range.front;
            /// It's presumed that chunks are sorted
            assert(next >= _front);
            /// Check if the chunk can be joined with the previous one
            if (_front.end >= next.beg) {
                /// update end of _front
                _front.end = max(_front.end, next.end);
                _range.popFront(); /// next is consumed
            } else {
                /// can't join
                break;
            }
        }
    }

    Chunk _front;
    bool _empty = false;
}

/// Params: 
///      r - range of chunks, sorted by leftmost coordinate
/// Returns:
///      range of non-overlapping chunks, covering the same subset
///      as original chunks
auto nonOverlappingChunks(R)(R r) 
    if (is(ElementType!R == Chunk)) 
{
    return NonOverlappingChunks!R(r);
}
