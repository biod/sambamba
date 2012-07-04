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

    int opCmp(VirtualOffset other) const nothrow {
        if (this.voffset > other.voffset) { return  1; }
        if (this.voffset < other.voffset) { return -1; }
        return 0;
    }

    int opCmp(ulong voffset) const nothrow {
        return opCmp(VirtualOffset(voffset));
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
}

