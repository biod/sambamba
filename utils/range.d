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
module utils.range;

import std.range;
import std.exception;
import std.algorithm;
import std.parallelism;
import std.functional;
import std.array;

/// Keeps a cyclic buffer of size $(D amount)
/// which is filled at the construction.
/// After that, each popFront() is accompanied
/// by fetching next element from the original range.
///
/// The function is useful when, for instance, range of Tasks
/// is being decorated, because it allows to keep a certain amount
/// of them being executed simultaneously, utilizing all
/// CPU cores.
auto prefetch(Range)(Range r, size_t amount) {

    enforce(amount > 0, "Amount of elements to prefetch must be positive");

    struct Result {
        alias ElementType!Range E;

        this(Range range, size_t amount) {
            _elements = new E[amount];
            _range = range;
            _amount = amount;
            foreach (i; 0 .. _amount) {
                if (_range.empty) {
                    break;
                }
                _elements[i] = _range.front;
                ++ _read;
                _range.popFront();
            }
        }
        
        bool empty() @property {
            return _range.empty() && _read == _consumed;
        }
        
        E front() @property {
            return _elements[_consumed % _amount];
        }

        void popFront() @property {
            if (_range.empty) {
                ++ _consumed;
                return;
            }

            _elements[_consumed % _amount] = _range.front;
            ++ _consumed;

            _range.popFront();
            ++ _read;
        }
    private:
        Range _range;
        size_t _read = 0;
        size_t _consumed = 0;
        size_t _amount;
        E[] _elements = void;
    }

    return Result(r, amount);
}

unittest {
    import std.algorithm;

    ubyte[] emptyrange = [];
    assert(equal(emptyrange, prefetch(emptyrange, 42)));

    auto range = [1, 2, 3, 4, 5];
    assert(equal(range, prefetch(range, 1)));
    assert(equal(range, prefetch(range, 3)));
    assert(equal(range, prefetch(range, 5)));
    assert(equal(range, prefetch(range, 7)));
}

/// Takes arbitrary input range as an input and returns
/// another range which produces arrays of original elements
/// of size $(D chunk_size).
///
/// Useful for setting granularity in parallel applications.
/// $(D std.algorithm.joiner) composed with $(D chunked) 
/// produces same elements as were in the original range.
/// 
/// The difference from $(D std.range.chunks) is that
/// any input range is allowed, no slicing or length is required.
/// The cost is memory allocations for chunks.
auto chunked(R)(R range, uint chunk_size) {

    alias ElementType!R E;

    struct Result {

        this(R range, uint chunk_size) {
            enforce(chunk_size > 0);
            this.range = range;
            this.chunk_size = chunk_size; 
            fillBuffer();
        }

        bool empty() @property {
            return buffer.length == 0;
        }

        E[] front() @property {
            return buffer;    
        }

        void popFront() {
            fillBuffer();
        }

    private:
        R range;
        uint chunk_size;

        E[] buffer;

        void fillBuffer() {
            buffer = uninitializedArray!(E[])(chunk_size);
            for (auto i = 0; i < chunk_size; i++) {
                if (range.empty) {
                    buffer.length = i;
                    break;
                }
                buffer[i] = range.front;
                range.popFront();
            }
        }
    }

    return Result(range, chunk_size);
}

unittest {
    import std.algorithm;

    assert(equal(chunked(iota(1, 6), 2), [[1, 2], [3, 4], [5]]));
    assert(equal(chunked(iota(1, 7), 2), [[1, 2], [3, 4], [5, 6]]));
    assert(equal(chunked([1], 10), [[1]]));
    assert(equal(chunked(iota(1, 10), 7), [[1, 2, 3, 4, 5, 6, 7], [8,9]]));

    auto r = iota(25);
    assert(equal(joiner(chunked(r, 7)), r));
} 

/// Version of parallel map using cyclic buffer with prefetching.
/// Uses combination of chunked, prefetch, joiner, and std.parallelism.
///
/// The analogue in Haskell is Control.Parallel.Strategies.parBuffer
/// 
/// Params:
///     prefetch_amount -   how many chunks will be prefetched
///     chunk_size      -   the maximum size of each chunk
auto parallelTransform(alias func, Range)(Range r, 
                                          uint chunk_size=1, 
                                          uint prefetch_amount=totalCPUs-1)
{
    alias ElementType!Range E;

    static auto createTask(E[] elements) {
        auto task =  task!(pipe!(map!(unaryFun!func), array))(elements);
        taskPool.put(task);
        return task;
    }

    if (prefetch_amount == 0) {
        prefetch_amount = 1;
    }

    auto chunks = chunked(r, chunk_size);
    auto tasks = map!createTask(chunks);
    auto prefetched = prefetch(tasks, prefetch_amount);
    return joiner(map!"a.yieldForce()"(prefetched));
}

unittest {
    auto range = iota(100);
    assert(equal(parallelTransform!"a * a"(range), map!"a * a"(range)));
} 
