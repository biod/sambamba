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
    }

    this(cram_fd* fd, bool reuse_buf, IBamSamReader reader,
				 cram_range* range=null) {
        _fd = fd;
				if (range !is null)
						_range = *range;
				else
						_range.refid = -4;
        _reader = reader;
        _reuse_buf = reuse_buf;
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
        auto cr = cram_get_seq(_fd);
        if (cr is null) {
            if (_fd.err == 0)
                _empty = true;
            else
                throw new CramException("Failed to read CRAM record");
        }
				if (_empty) return;
				_front = cram2bam(_fd.header, _fd, _fd.ctr.slice, cr,
						_fd.ctr.curr_rec - 1);
        _front.associateWithReader(_reader);
    }

		char[1024] name_buf;

    private BamRead cram2bam(SAM_hdr *bfd, cram_fd *fd, cram_slice *s,
            cram_record *cr, int rec) {
        char* name;
        int name_len;
        if (cr.name_len > 0) {
            name = cast(char*)(s.name_blk.data + cr.name);
            name_len = cr.name_len;
        } else {
            auto p = name_buf.ptr;
            name = p;
            auto prefix = fd.prefix;
            while (*prefix != '\0')
                *p++ = *prefix++;
            bio.core.utils.format.write(p, ':');
            bio.core.utils.format.write(p, s.id);
            bio.core.utils.format.write(p, ':');
            if (cr.mate_line >= 0 && cr.mate_line < rec)
                bio.core.utils.format.write(p, cr.mate_line);
            else
                bio.core.utils.format.write(p, rec);
            name_len = cast(int)(p - name);
        }

        if (cr.rg < -1 || cr.rg >= bfd.nrg)
            throw new CramException("Read group id is out of range");

        // extra 4 bytes account for RGZ and zero
        auto rg_len = (cr.rg != -1) ? (bfd.rg[cr.rg].name_len + 4) : 0;

        enforce(s.seqs_blk.data !is null && s.qual_blk.data !is null);

        auto name_offset = 32;
        auto cigar_offset = name_offset + name_len + 1;
        auto seq_offset = cigar_offset + cr.ncigar * uint.sizeof;
        auto qual_offset = seq_offset + (cr.len + 1) / 2;
        auto tag_offset = qual_offset + cr.len;

        ubyte[] raw_data = allocate(cr.aux_size + rg_len + tag_offset);
        uint bin_mq_nl = (reg2bin(cr.apos - 1, cr.aend) << 16) |
            (cr.mqual << 8) | (name_len + 1);
        uint flag_nc = (cr.flags << 16) | cr.ncigar;

        *cast(int*)(raw_data.ptr) = cr.ref_id;
        *cast(int*)(raw_data.ptr + 4) = cr.apos - 1;
        *cast(uint*)(raw_data.ptr + 8) = bin_mq_nl;
        *cast(uint*)(raw_data.ptr + 12) = flag_nc;
        *cast(uint*)(raw_data.ptr + 16) = cast(uint)cr.len;
        *cast(int*)(raw_data.ptr + 20) = cr.mate_ref_id;
        *cast(int*)(raw_data.ptr + 24) = cr.mate_pos - 1;
        *cast(int*)(raw_data.ptr + 28) = cr.tlen;

        raw_data[name_offset .. name_offset+name_len] = (cast(ubyte[])name[0 .. name_len])[];
        auto cigar_ptr = cast(uint*)(raw_data.ptr + cigar_offset);
        cigar_ptr[0 .. cr.ncigar] = (s.cigar + cr.cigar)[0 .. cr.ncigar];
        auto src = cast(string)((s.seqs_blk.data + cr.seq)[0 .. cr.len]);
        auto dst = raw_data.ptr + seq_offset;
        foreach (k; 0 .. cr.len / 2) {
            auto b1 = Base16(src[2*k]).internal_code;
            auto b2 = Base16(src[2*k+1]).internal_code;
            dst[k] = cast(ubyte)((b1 << 4) | b2);
        }
        if ((cr.len & 1) != 0)
            dst[cr.len / 2] = cast(ubyte)(Base16(src[$-1]).internal_code << 4);

        auto quals = (s.qual_blk.data + cr.qual)[0 .. cr.len];
        raw_data[qual_offset .. $][0 .. cr.len] = quals[];

        auto tag_data = raw_data[tag_offset .. $][0 .. cr.aux_size];
        tag_data[] = (s.aux_blk.data + cr.aux)[0 .. cr.aux_size];

        if (cr.rg != -1) {
            auto rg_name_str = bfd.rg[cr.rg].name[0 .. rg_len - 4];
            auto rg_name = Value(cast(string)rg_name_str);
            (tag_data.ptr + cr.aux_size)[0 .. 2] = cast(ubyte[])"RG";
            emplaceValue(tag_data.ptr + cr.aux_size + 2, rg_name);
        }

        return BamRead(raw_data);
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

