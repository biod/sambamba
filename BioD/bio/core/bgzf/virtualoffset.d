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
module bio.core.bgzf.virtualoffset;

import std.conv;

/// Structure representing virtual offset in BGZF-compressed file.
struct VirtualOffset {
    /// Params:
    ///     
    ///     coffset =    unsigned byte offset into the BGZF file 
    ///                  to the beginning of a BGZF block.
    ///                  Must be strictly less than 2^48.
    ///
    ///     uoffset =    unsigned byte offset into the uncompressed
    ///                  data stream represented by that BGZF block
    this(ulong coffset, ushort uoffset) nothrow @safe
    in {
        assert(coffset < (1UL<<48));
    } 
    body {
        voffset = (coffset << 16) | uoffset;
    }
    
    /// Set both coffset and uoffset packed as (coffset<<16)|uoffset
    this(ulong voffset) nothrow @safe {
        this.voffset = voffset;
    }

    /// ditto
    ulong coffset() @property const nothrow @safe pure {
        return voffset >> 16;
    }
    
    /// ditto
    ushort uoffset() @property const nothrow @safe pure {
        return voffset & 0xFFFF;
    }

    int opCmp(const ref VirtualOffset other) const nothrow @safe pure {
        if (this.voffset > other.voffset) { return  1; }
        if (this.voffset < other.voffset) { return -1; }
        return 0;
    }

    bool opEquals(const ref VirtualOffset other) const nothrow @safe {
        return this.voffset == other.voffset;
    }

    bool opEquals(ulong voffset) const nothrow @safe {
        auto vo = VirtualOffset(voffset);
        return opEquals(vo);
    }

    ulong opCast() const nothrow @safe pure {
        return voffset;
    }

    /// String representation in format "<coffset>/<uoffset>"
    string toString() {
        return to!string(coffset) ~ "/" ~ to!string(uoffset);
    }

private:
    ulong voffset;
}

unittest {
    auto voffset = VirtualOffset(100500, 42);
    assert(voffset.coffset == 100500);
    assert(voffset.uoffset == 42);
    assert(voffset == (100500UL << 16) + 42UL);
    assert(cast(ulong)voffset == (100500UL << 16) + 42UL);
}
