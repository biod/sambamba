module cram.slicereader;

import cram.htslib, cram.wrappers, cram.exception;
import bio.bam.abstractreader, bio.bam.read;

alias AllocateFunc = ubyte[] function(size_t);

ubyte[] reusableChunkAllocator(size_t size) {
    static ubyte[] _block;
    if (size > _block.length) {
        _block.assumeSafeAppend();
        _block.length = size;
    }
    return _block[0 .. size];
}

ubyte[] chunkAllocator(size_t size) {
    static ubyte[] _block;
    static size_t _capa = 65536;
    static size_t _used;
            
    if (size + _used <= _capa && _block !is null) {
        auto result = _block[_used .. $][0 .. size];
        _used += size;
        return result;
    } else {
        if (size > _capa)
            _capa = size;
        _block = uninitializedArray!(ubyte[])(_capa);
        _used = size;
        return _block[0 .. size];
    }
}

AllocateFunc bamReadAlloc(bool buffer_is_safe_to_reuse) {
    if (buffer_is_safe_to_reuse)
        return &reusableChunkAllocator;
    else
        return &chunkAllocator;
}

struct CramSliceReader {
    private {
        CramSlice _slice;
        int _n_curr;
        size_t _n_total;
        alias RcBamSeq = RcPtr!(bam_seq_t, bam_destroy1);
        RcBamSeq _bam_seq_ptr;
        IBamSamReader _reader;
        AllocateFunc _alloc;
    }

    BamRead front;
    bool empty;

    this(CramSlice slice, IBamSamReader reader, AllocateFunc alloc) {
        _bam_seq_ptr = RcBamSeq(bam_init1());
        _slice = slice;
        _reader = reader;
        _alloc = alloc;
        _n_total = slice.hdr.num_records;
        popFront();
    }

    void popFront() {
        if (_n_curr == _n_total) {
            empty = true;
            return;
        }

        int ret = cram_to_bam(_slice.fd.header,
                              _slice.fd.data.ptr,
                              _slice.data.ptr,
                              &_slice.crecs[_n_curr],
                              _n_curr,
                              &_bam_seq_ptr.data.ptr);
        if (ret == -1)
            throw new CramException("Failure in cram_to_bam");

        front = toBamRead(_bam_seq_ptr);
        front.associateWithReader(_reader);

        ++_n_curr;
    }

    private BamRead toBamRead(bam_seq_t* record) {
        auto chunk = _alloc(32 + record.l_data);
        chunk[32 .. $] = record.data[0 .. record.l_data];
        auto ptr = cast(void*)(chunk.ptr);
        *(cast(int*)ptr + 0) = record.core.tid;
        *(cast(int*)ptr + 1) = record.core.pos;

        auto bin = record.core.bin;
        auto mq = record.core.qual;
        auto nl = record.core.l_qname;
        *(cast(uint*)ptr + 2) = (bin << 16) | (mq << 8) | nl;
        auto flag = record.core.flag;
        auto nc = record.core.n_cigar;
        *(cast(uint*)ptr + 3) = (flag << 16) | nc;

        *(cast(int*)ptr + 4) = record.core.l_qseq;
        *(cast(int*)ptr + 5) = record.core.mtid;
        *(cast(int*)ptr + 6) = record.core.mpos;
        *(cast(int*)ptr + 7) = record.core.isize;
        return BamRead(chunk);
    }
}

auto bamReads(CramSlice slice, IBamSamReader reader=null, 
              AllocateFunc alloc=&chunkAllocator)
{
    return CramSliceReader(slice, reader, alloc); 
}
