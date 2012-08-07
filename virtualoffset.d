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
module virtualoffset;

import std.conv;

/// Structure representing virtual offset in BAM file.
struct VirtualOffset {
    /// Params:
    ///     
    ///     coffset =    unsigned byte offset into the BGZF file 
    ///                  to the beginning of a BGZF block.
    ///                  Must be strictly less than 2^48.
    ///
    ///     uoffset =    unsigned byte offset into the uncompressed
    ///                  data stream represented by that BGZF block
    this(ulong coffset, ushort uoffset) nothrow 
    in {
        assert(coffset < (1UL<<48));
    } 
    body {
        voffset = (coffset << 16) | uoffset;
    }
    
    /// Set both coffset and uoffset packed as (coffset<<16)|uoffset
    this(ulong voffset) nothrow {
        this.voffset = voffset;
    }

    /// ditto
    ulong coffset() @property const nothrow {
        return voffset >> 16;
    }
    
    /// ditto
    ushort uoffset() @property const nothrow {
        return voffset & 0xFFFF;
    }

    int opCmp(const ref VirtualOffset other) const nothrow {
        if (this.voffset > other.voffset) { return  1; }
        if (this.voffset < other.voffset) { return -1; }
        return 0;
    }

    bool opEquals(const ref VirtualOffset other) const nothrow {
        return this.voffset == other.voffset;
    }

    bool opEquals(ulong voffset) const nothrow {
        return opEquals(VirtualOffset(voffset));
    }

    ulong opCast() const nothrow {
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
