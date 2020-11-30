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
module bio.std.hts.utils.array;

import core.stdc.string;
import std.traits;

/// Modifies array in-place so that $(D slice) is replaced by
/// $(D replacement[]).
///
/// WARNING: it's you who is responsible that $(D slice) is indeed
/// a slice of $(D s).
void replaceSlice(T, U)(ref T[] s, in U[] slice, in T[] replacement)
    if (is(Unqual!U == T)) 
{

    auto offset = slice.ptr - s.ptr;
    auto slicelen = slice.length;
    auto replen = replacement.length;

    auto newlen = s.length - slicelen + replen;

    if (slicelen == replen) {
        s[offset .. offset + slicelen] = replacement;
        return;
    }

    if (replen < slicelen) {
        // overwrite piece of slice
        s[offset .. offset + replen] = replacement;
        // and then move the remainder
        memmove(s.ptr + (offset + replen),
                s.ptr + (offset + slicelen),
                (newlen - offset - replen) * T.sizeof);

        s.length = newlen;
        return;
    }

    // replen > slicelen
    s.length = newlen;
    // here, first move the remainder
    memmove(s.ptr + (offset + replen),
            s.ptr + (offset + slicelen),
            (newlen - offset - replen) * T.sizeof);
    // and then overwrite
    s[offset .. offset + replen] = replacement;
}

/// Does almost the same, but does not require $(D replacement),
/// instead only its length, $(D n) bytes. This is useful for
/// avoiding memory allocations.
void prepareSlice(T, U)(ref T[] s, in U[] slice, size_t n)
    if (is(Unqual!U == T))
{

    auto offset = slice.ptr - s.ptr;
    auto slicelen = slice.length;
    auto replen = n;

    auto newlen = s.length - slicelen + replen;

    if (slicelen == replen) {
        return;
    }

    if (replen < slicelen) {
        memmove(s.ptr + (offset + replen),
                s.ptr + (offset + slicelen),
                (newlen - offset - replen) * T.sizeof);

        s.length = newlen;
        return;
    }

    // replen > slicelen
    s.length = newlen;
    memmove(s.ptr + (offset + replen),
            s.ptr + (offset + slicelen),
            (newlen - offset - replen) * T.sizeof);
}


unittest {
    auto array = [1, 2, 3, 4, 5];
    replaceSlice(array, array[2 .. 4], [-1, -2, -5]);
    assert(array == [1, 2, -1, -2, -5, 5]);
    replaceSlice(array, array[1 .. 4], cast(int[])[]);
    assert(array == [1, -5, 5]);
    replaceSlice(array, array[0 .. 1], [3]);
    assert(array == [3, -5, 5]);
}
