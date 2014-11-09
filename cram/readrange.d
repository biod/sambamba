module cram.readrange;

import cram.htslib;
import cram.exception;

import bio.bam.abstractreader, bio.bam.read, bio.core.utils.format,
			 bio.core.base, bio.bam.utils.value, bio.bam.tagvalue, bio.core.base,
       bio.bam.bai.bin;
import std.array;

struct CramBamReadRange {
    private {
        cram_fd* _fd;
        cram_range _range;
        IBamSamReader _reader;
        bool _popped_first;

        bam_seq_t* _bam_seq_ptr;
    }

    this(cram_fd* fd, bool reuse_buf, IBamSamReader reader,
         cram_range* range=null)
    {
        _fd = fd;
        if (range !is null)
            _range = *range;
        else
            _range.refid = -4;
        _reader = reader;
        _reuse_buf = reuse_buf;
        _bam_seq_ptr = bam_init1();
    }

    this(this) {
        _bam_seq_ptr = bam_init1();
    }

    ~this() {
        bam_destroy1(_bam_seq_ptr);
    }

    private void popFrontIfNeeded() {
        if (!_popped_first) {
            if (_range.refid != -4) {
                cram_set_option(_fd, cram_option.CRAM_OPT_RANGE, &_range);
            }
            popFront();
            _popped_first = true;
        }
    }

    private bio.bam.read.BamRead _front;
    bio.bam.read.BamRead front() @property {
        popFrontIfNeeded();
        return _front;
    }
    private bool _empty;
    bool empty() @property {
        popFrontIfNeeded();
        return _empty;
    }

    void popFront() {
        int ret = cram_get_bam_seq(_fd, &_bam_seq_ptr);
        if (ret == -1) {
            if (_fd.err == 0)
                _empty = true;
            else
                throw new CramException("Failed to read CRAM record");
        }
        if (_empty) return;
        _front = toBamRead(_bam_seq_ptr);
        _front.associateWithReader(_reader);
    }

    char[1024] name_buf;

    private BamRead toBamRead(bam_seq_t* record) {
        auto chunk = allocate(32 + record.l_data);
        chunk[32 .. $] = record.data[0 .. record.l_data];
        auto ptr = cast(void*)(chunk.ptr);
        *(cast(int*)ptr + 0) = record.core.tid;
        *(cast(int*)ptr + 1) = record.core.pos;

        // UNDEFINED BEHAVIOUR (bitfields)
        auto bin = record.core.bin_mq_nl & 0xFFFF;
        auto mq = (record.core.bin_mq_nl >> 16) & 0xFF;
        auto nl = record.core.bin_mq_nl >> 24;
        *(cast(uint*)ptr + 2) = (bin << 16) | (mq << 8) | nl;
        auto flag = record.core.flag_nc & 0xFFFF;
        auto nc = record.core.flag_nc >> 16;
        *(cast(uint*)ptr + 3) = (flag << 16) | nc;

        *(cast(int*)ptr + 4) = record.core.l_qseq;
        *(cast(int*)ptr + 5) = record.core.mtid;
        *(cast(int*)ptr + 6) = record.core.mpos;
        *(cast(int*)ptr + 7) = record.core.isize;
        return BamRead(chunk);
    }

    private {
        bool _reuse_buf;
        size_t _used;
        ubyte[] _block;
        size_t _capa = 65536;
        ubyte[] allocate(size_t size) {
            if (_reuse_buf) {
                if (size > _block.length) {
                  _block.assumeSafeAppend();
                  _block.length = size;
                }
                return _block[0 .. size];
            }

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
    }
}

