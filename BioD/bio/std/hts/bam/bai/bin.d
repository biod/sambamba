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
module bio.std.hts.bam.bai.bin;

import bio.core.bgzf.chunk;
import bio.std.hts.bam.constants;

/// Distinct bin
struct Bin {

    /// Construct a bin with an id
    this(uint id) nothrow {
        this.id = id;
    }

    uint id; /// bin number
    Chunk[] chunks; 

    /// How deep the bin is in the tree
    int level() @property const nothrow {
        if (id == 0) return 0;
        if (id < 9) return 1;
        if (id < 73) return 2;
        if (id < 585) return 3;
        if (id < 4681) return 4;
        return 5;
    }

    /// Returns whether the bin is a leaf in the B-tree
    bool is_leaf() @property const nothrow {
        return id > BAI_MAX_NONLEAF_BIN_ID;
    }

    /// Check if bin can overlap with a region
    bool canOverlapWith(int begin, int end) const nothrow {
        if (id == 0) return true;
        if (id > BAI_MAX_BIN_ID) return false;

        /// The following code is based on reg2bins() function
        if (begin < 0) begin = 0;
        auto magic_number = 4681;
        auto b = begin >> 14;
        auto e = end   >> 14;

        while (true) {
            auto delta = id - magic_number;
            if (b <= delta && delta <= e) return true;

            magic_number >>= 3;

            if (magic_number == 0) return false;

            b >>= 3;
            e >>= 3;
        }
    } 
}

/// Returns bin number for [beg, end) interval (zero-based).
/// Taken from SAM/BAM specification.
ushort reg2bin(int beg, int end) {
    if (end == beg) end = beg + 1; // edge case

    --end;
    if (beg>>14 == end>>14) return cast(ushort)(((1<<15)-1)/7 + (beg>>14));
    if (beg>>17 == end>>17) return cast(ushort)(((1<<12)-1)/7 + (beg>>17));
    if (beg>>20 == end>>20) return cast(ushort)(((1<<9)-1)/7  + (beg>>20));
    if (beg>>23 == end>>23) return cast(ushort)(((1<<6)-1)/7  + (beg>>23));
    if (beg>>26 == end>>26) return cast(ushort)(((1<<3)-1)/7  + (beg>>26));
    return 0;
}
