module bai.bin;

import bai.chunk;

/// Distinct bin
struct Bin {
    uint id; /// bin number
    Chunk[] chunks; 

    /// How deep the bin is in the tree
    int level() @property {
        if (id == 0) return 0;
        if (id < 9) return 1;
        if (id < 73) return 2;
        if (id < 585) return 3;
        if (id < 4681) return 4;
        return 5;
    }

	/// Check if bin can overlap with a region
	bool canOverlapWith(int begin, int end) {
		if (id == 0) return true;

		/// The following code is based on reg2bins() function
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
