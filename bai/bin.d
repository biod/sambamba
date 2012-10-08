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
module bai.bin;

import bai.chunk;
import constants;

/// Distinct bin
struct Bin {

    /// Construct a bin with an id
    this(ushort id) nothrow {
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
    --end;
    if (beg>>14 == end>>14) return cast(ushort)(((1<<15)-1)/7 + (beg>>14));
    if (beg>>17 == end>>17) return cast(ushort)(((1<<12)-1)/7 + (beg>>17));
    if (beg>>20 == end>>20) return cast(ushort)(((1<<9)-1)/7  + (beg>>20));
    if (beg>>23 == end>>23) return cast(ushort)(((1<<6)-1)/7  + (beg>>23));
    if (beg>>26 == end>>26) return cast(ushort)(((1<<3)-1)/7  + (beg>>26));
    return 0;
}
