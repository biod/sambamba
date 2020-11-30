/*
    This file is part of BioD.
    Copyright (C) 2013-2017    Artem Tarasov <lomereiter@gmail.com>

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
module bio.core.utils.outbuffer;

import std.array;
import std.exception;

///
class OutBuffer {
    private {
        ubyte[] _heap;

        ubyte* _heap_ptr() @property { return _heap.ptr; }
        size_t _heap_used;
        size_t _heap_capacity;
    }

    ///
    this(size_t initial_capacity) {
        _heap = uninitializedArray!(ubyte[])(initial_capacity);
        _heap_capacity = initial_capacity;
    }

    ///
    ubyte[] data() @property {
        return _heap_ptr[0 .. _heap_used];
    }

    ///
    size_t length() @property const {
        return _heap_used;
    }

    /// Remove last elements such that new size is equal to $(D size).
    void shrink(size_t size) {
        enforce(size <= length);
        _heap_used = size;
    }

    ///
    size_t capacity() @property const {
        return _heap_capacity;
    }
    
    /// ditto
    void capacity(size_t new_capacity) @property {
        if (new_capacity <= _heap_capacity)
            return;

        _heap.length = new_capacity;
        _heap_capacity = new_capacity;
    }
        
    ///
    void put(T)(T bytes) if (is(T == ubyte[])) {
        size_t needed = bytes.length + _heap_used;
        if (needed > _heap_capacity) {
            do {
                _heap_capacity = _heap_capacity * 3 / 2;
            } while (_heap_capacity < needed);
            _heap.length = _heap_capacity;
        }
        _heap_ptr[_heap_used .. _heap_used + bytes.length] = bytes[];
        _heap_used = needed;
    }

    /// Dumps raw bytes into the buffer. No endianness conversion or whatever.
    void put(T)(auto ref T value) if (!is(T == ubyte[])) {
        put((cast(ubyte*)(&value))[0 .. T.sizeof]);
    }

    /// Responsibility that there's enough capacity is on the user
    void putUnsafe(T)(T bytes) if (is(T == ubyte[])) {
        _heap_ptr[_heap_used .. _heap_used + bytes.length] = bytes[];
        _heap_used += bytes.length;
    }

    /// ditto
    void putUnsafe(T)(auto ref T value) if (!is(T == ubyte[])) {
        putUnsafe((cast(ubyte*)(&value))[0 .. T.sizeof]);
    }

    ///
    void clear() {
        _heap_used = 0;
    }
}
