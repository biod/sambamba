/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.core.utils.range;

import bio.core.utils.roundbuf;

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
            _roundbuf = RoundBuf!E(amount);
            _range = range;
            foreach (i; 0 .. amount) {
                if (_range.empty) {
                    break;
                }
                _roundbuf.put(_range.front);
                _range.popFront();
            }
        }
        
        bool empty() @property {
            return _range.empty && _roundbuf.empty;
        }
        
        auto ref front() @property {
            return _roundbuf.front;
        }

        void popFront() @property {
            assert(!_roundbuf.empty);

            if (_range.empty) {
                _roundbuf.popFront();
                return;
            }

            _roundbuf.popFront();
            _roundbuf.put(_range.front);

            _range.popFront();
        }
    private:
        Range _range;
        RoundBuf!E _roundbuf;
    }

    return Result(r, amount);
}

///
struct Cached(R) {
    private {
        alias ElementType!R E;
        R _range;
        E _front;
        bool _empty;
    }

    this(R range) {
        _range = range;
        popFront();
    }

    auto front() { return _front; }
    bool empty() { return _empty; }
    void popFront() {
        if (_range.empty) {
            _empty = true;
        } else {
            _front = _range.front;
            _range.popFront();
        }
    }
}

/// Caches front element.
auto cached(R)(R range) {
    return Cached!R(range);
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

    assert(equal(range, cached(range)));
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

struct PrefixSum(S) {
    private {
        S _sequence;
        ElementType!S _sum;
    }

    this(S sequence) {
        _sequence = sequence;
        if (!_sequence.empty) {
            _sum = _sequence.front;
        }
    }

    bool empty() @property {
        return _sequence.empty;
    }

    ElementType!S front() @property {
        return _sum;
    }

    void popFront() {
        _sequence.popFront();
        if (!_sequence.empty) {
            _sum += _sequence.front;
        }
    }
}

/// Prefix sum.
PrefixSum!S prefixSum(S)(S sequence) {
    return PrefixSum!S(sequence);
}

unittest {
    auto range = iota(5);
    assert(equal(prefixSum(range), [0, 1, 3, 6, 10]));
}
