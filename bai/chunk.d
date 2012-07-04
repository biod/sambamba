module bai.chunk;

import virtualoffset;

struct Chunk {
    VirtualOffset beg; /// virtual file offset of the start of the chunk
    VirtualOffset end; /// virtual file offset of the end of the chunk

	/// Compare beginnings
	int opCmp(Chunk other) const nothrow {
		if (beg < other.beg) return -1;
		if (beg > other.beg) return 1;
		if (end < other.end) return -1;
		if (end > other.end) return 1;
		return 0;
	}
}
