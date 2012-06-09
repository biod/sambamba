module utils.array;

import std.c.string;

/// Modifies array in-place so that $(D slice) is replaced by
/// $(D replacement[]).
///
/// WARNING: it's you who is responsible that $(D slice) is indeed
/// a slice of $(D s).
void replaceSlice(T)(ref T[] s, in T[] slice, in T[] replacement)
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

unittest {
    auto array = [1, 2, 3, 4, 5];
    replaceSlice(array, array[2 .. 4], [-1, -2, -5]);
    assert(array == [1, 2, -1, -2, -5, 5]);
    replaceSlice(array, array[1 .. 4], cast(int[])[]);
    assert(array == [1, -5, 5]);
    replaceSlice(array, array[0 .. 1], [3]);
    assert(array == [3, -5, 5]);
}
