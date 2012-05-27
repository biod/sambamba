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
