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
module bio.core.utils.roundbuf;

import std.exception;

/// Cyclic buffer, see https://en.wikipedia.org/wiki/Circular_buffer
/// Bails out if you try to push beyond its size
///
/// Note the RoundBuf has copy semantics. So if you pass in a
/// Struct it will copy its contents.

struct RoundBuf(T) {

    private {
      T[] _items = void; // max_length is the size of round buf
      size_t _put;       // moves the counter forward
      size_t _start;     // current start of the buffer
    }

    /** initializes round buffer of size $(D n) */
    this(size_t n) {
        _items = new T[n];
    }

    /// Input range primitives
    bool empty() @property const {
        return _put == _start;
    }

    /// Get the item at the front (non-destructive)
    auto ref front() @property {
        enforce(!empty, "roundbuffer is empty");
        return _items[_start % $];
    }

    /// Move the front forward (destructive)
    void popFront() {
        enforce(!empty, "roundbuffer is empty");
        ++_start;
    }

    /// Returns the tail item (non-destructive)
    auto ref back() @property {
      enforce(!empty, "roundbuffer is empty");
      return _items[(_put - 1) % $];
    }

    /// Add and item at the tail and move pointer forward (destructive)
    void put(T item) {
        enforce(!full, "roundbuffer is full");
        enforce(_put < _put.max, "ringbuffer size_t overflow");
        _items[_put % $] = item;
        ++_put;
    }

    /// Check if buffer is full
    bool full() @property const {
      return _put == _start + max_length;
    }

    /// Current number of elements
    size_t length() @property const {
        return _put - _start;
    }

    size_t max_length() @property const {
      return _items.length;
    }
}

unittest {
  import std.stdio;
    auto buf = RoundBuf!int(4);
    assert(buf.empty);

    buf.put(1);
    buf.put(2);
    assert(buf.length == 2);
    assert(buf.front == 1);
    buf.popFront();
    buf.put(1);
    buf.put(0);
    buf.put(3);
    assert(buf.full);
    buf.popFront();
    buf.put(4);
    buf.popFront();
    buf.popFront();
    assert(buf.front == 3);
    buf.popFront();
    assert(buf.front == 4);
    buf.put(4);
    buf.put(5);
    buf.put(6);
    assert(buf.length == buf.max_length);
    // should bomb out
    assertThrown!Exception(buf.put(7));
}
