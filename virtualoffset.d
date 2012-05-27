module virtualoffset;

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
    this(ulong coffset, ushort uoffset) 
    in {
        assert(coffset < (1UL<<48));
    } 
    body {
        voffset = (coffset << 16) | uoffset;
    }
    
    /// Set both coffset and uoffset packed as (coffset<<16)|uoffset
    this(ulong voffset) {
        this.voffset = voffset;
    }

    /// ditto
    ulong coffset() @property {
        return voffset >> 16;
    }
    
    /// ditto
    ushort uoffset() @property {
        return voffset & 0xFFFF;
    }

    int opCmp(VirtualOffset other) {
        if (this.voffset > other.voffset) { return  1; }
        if (this.voffset < other.voffset) { return -1; }
        return 0;
    }

    int opCmp(ulong voffset) {
        return opCmp(VirtualOffset(voffset));
    }

private:
    ulong voffset;
}

unittest {
    auto voffset = VirtualOffset(100500, 42);
    assert(voffset.coffset == 100500);
    assert(voffset.uoffset == 42);
}

