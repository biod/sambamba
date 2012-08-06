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
module utils.array;

import std.c.string;
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
