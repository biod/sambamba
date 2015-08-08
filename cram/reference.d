module cram.reference;

import cram.exception, cram.wrappers, cram.slicereader;
import cram.htslib;

import bio.bam.abstractreader;
import bio.bam.referenceinfo;
import bio.bam.read;
import std.string;
import std.parallelism;
import std.range, std.algorithm;

int ref_seq_id(BamRead r) { return r.ref_id; }
int ref_seq_start(BamRead r) { return r.position + 1; }
int ref_seq_span(BamRead r) { return r.basesCovered(); }

struct PositionChecker {
    cram_range range;

    private CramFilterResult passImpl(T)(T x) {
        if (x.ref_seq_id < range.refid) return CramFilterResult.skip;
        if (x.ref_seq_id > range.refid) return CramFilterResult.stop;
        if (x.ref_seq_start > range.end) return CramFilterResult.stop;

        // order of checks matters: ref_seq_span is relatively slow for reads
        if (x.ref_seq_start + x.ref_seq_span - 1 < range.start) 
            return CramFilterResult.skip;
        
        return CramFilterResult.pass;
    }

    auto pass(cram_container* c) {
        if (c.length == 0) return CramFilterResult.skip;
        return passImpl(c);
    }

    auto pass(cram_slice* s) { return passImpl(s.hdr); }
    auto pass(BamRead r) { return passImpl(r); }
}

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

    private TaskPool _task_pool;

    this(IBamSamReader reader, CramFd fd, bool seq_op, 
         int ref_id, ReferenceSequenceInfo info,
         TaskPool task_pool)
    {
        _reader = reader;
        _fd = fd;
        _seq_op = seq_op;
        _ref_id = ref_id;
        _info = info;
        _task_pool = task_pool;
    }

    /// Get alignments overlapping [start, end) region.
    /// $(BR)
    /// Coordinates are 0-based.
    auto opSlice(uint start, uint end) {
        enforce(start < end, "start must be less than end");
        enforce(_ref_id >= 0, "invalid reference id");
        enforce(cram_index_load(_fd, toStringz(_reader.filename)) == 0,
                "couldn't load CRAM index");

        cram_range r;
        r.refid = _ref_id;
        r.start = start + 1;
        r.end = end;
        if (cram_seek_to_refpos(_fd, &r) == -1)
            throw new CramException("Failure in cram_seek_to_refpos");
        auto checker = PositionChecker(r);
        auto alloc = bamReadAlloc(_seq_op);

        return _fd.slices(c => checker.pass(c), s => checker.pass(s),
                          _task_pool)
                  .zip(repeat(checker), repeat(_reader), repeat(alloc))
                  .map!(x => bamReads(x[0], x[2], x[3])
                             .zip(repeat(x[1]))
                             .filter!(y => y[1].pass(y[0]) == CramFilterResult.pass)
                             .map!(y => y[0]))
                  .joiner2;
    }

    /// All alignments
    auto opSlice() {
        assert(length > 0);
        return opSlice(0, length);
    }

    private:
    IBamSamReader _reader;
    CramFd _fd;
    int _ref_id;
    bool _seq_op;
    ReferenceSequenceInfo _info;
}
