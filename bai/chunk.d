module bai.chunk;

struct Chunk {
    ulong beg; /// virtual file offset of the start of the chunk
    ulong end; /// virtual file offset of the end of the chunk

	/// Compare beginnings
	int opCmp(ref const Chunk other) const {
		if (beg < other.beg) return -1;
		if (beg > other.beg) return 1;
		if (end < other.end) return -1;
		if (end > other.end) return 1;
		return 0;
	}
}
