module cram.reference;

import cram.htslib, cram.exception, cram.readrange;
import bio.bam.abstractreader;
import bio.bam.referenceinfo;
import std.string;

struct ReferenceSequence {
    /// Name
    string name() @property const {
        return _info.name;
    }

    /// Length in base pairs
    int length() @property const {
        return _info.length;
    }

    /// Reference ID
    int id() @property const {
        return _ref_id;
    }

		this(IBamSamReader reader, cram_fd* fd, bool seq_op, 
				int ref_id, ReferenceSequenceInfo info) 
		{
				_reader = reader;
			  _fd = fd;
				_seq_op = seq_op;
				_ref_id = ref_id;
				_info = info;
		}

		/// Get alignments overlapping [start, end) region.
    /// $(BR)
    /// Coordinates are 0-based.
    auto opSlice(uint start, uint end) {
        enforce(start < end, "start must be less than end");
        enforce(_ref_id >= 0, "invalid reference id");
				import std.stdio;
				enforce(cram_index_load(_fd, toStringz(_reader.filename)) == 0,
						"couldn't load CRAM index");
				
				cram_range r;
				r.refid = _ref_id;
				r.start = start + 1;
				r.end = end;
				return CramBamReadRange(_fd, _seq_op, _reader, &r);
    }

		/// All alignments
		auto opSlice() {
			return opSlice(0, length);
		}

private:
		IBamSamReader _reader;
		cram_fd* _fd;
		int _ref_id;
		bool _seq_op;
		ReferenceSequenceInfo _info;
}
