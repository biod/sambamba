module cram.reader;

// uses some internal functions of htslib => static linking is necessary

import cram.htslib;
import cram.reference;
import cram.exception;
import cram.wrappers;
import cram.slicereader;

import bio.bam.abstractreader, bio.bam.referenceinfo,
       bio.sam.header, bio.bam.reference, bio.bam.read;
import std.string, std.array, std.range, std.algorithm, std.conv, std.typecons;
import std.parallelism;

class CramReader : IBamSamReader {
    private {
        string _fn;
        string _mode;
        CramFd _fd;
        TaskPool _task_pool;

        SamHeader _header; // (RE-)parsed lazily

        SamHeader parseHeader() {
            auto ptr = sam_hdr_str(_fd.header);
            auto len = sam_hdr_length(_fd.header);
            // FIXME: avoid reparsing header; turn it into an interface
            // with two implementations?
            return new SamHeader(cast(string)(ptr[0 .. len]));
        }

        ReferenceSequenceInfo[] _reference_sequences;
        int[string] _reference_sequence_dict;
    }

    this(string filename, TaskPool pool=taskPool) {
        _fn = filename;
        _task_pool = pool;

        _fd = openCram(filename);
        _reference_sequences.length = _fd.header.nref;
        foreach (k; 0 .. _fd.header.nref) {
            auto seq = _fd.header.ref_[k];
            _reference_sequences[k] = ReferenceSequenceInfo(seq.name.to!string,
                    seq.len);
            _reference_sequence_dict[_reference_sequences[k].name] = cast(int)k;
        }
    }

    ///
    bool hasReference(string reference) {
        return null != (reference in _reference_sequence_dict);
    }

    ///
    string filename() @property const {
        return _fn;
    }

    ///
    cram.reference.ReferenceSequence opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " is not present in the header");
        auto ref_id = _reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    ///
    cram.reference.ReferenceSequence reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return cram.reference.ReferenceSequence(this, fd(), _seq_op, ref_id,
                                                _reference_sequences[ref_id],
                                                _task_pool);
    }

    BamRead[] unmappedReads() {
        throw new Exception("* region unimplemented for CRAM files");
    }

    SamHeader header() @property {
        if (_header is null)
            _header = parseHeader();
        return _header;
    }

    const(ReferenceSequenceInfo)[] reference_sequences() @property const nothrow {
        return _reference_sequences;
    }

    private bool _seq_op;
    void assumeSequentialProcessing() {
        _seq_op = true;
    }

    private CramFd fd() {
        return _fn == "-" ? _fd : openCram(_fn);
    }

    auto reads() {
        alias R = CramFilterResult;
        auto s = fd().slices(c => c.length > 0 ? R.pass : R.skip,
                             null,
                             _task_pool);
        // no trust for delegates implementation => use zip-repeat trick
        auto reads = s.zip(repeat(this), repeat(bamReadAlloc(_seq_op)))
                      .map!(x => bamReads(x[0], x[1], x[2]))
                      .joiner2;
        return reads;
    }

    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads());
    }

    void createIndex(string fn_prefix=null) {
        int ret = cram_index_build(_fd, toStringz(fn_prefix is null ? _fn : fn_prefix));
        if (ret != 0) {
            throw new Exception("failed to build index for CRAM file " ~ _fn);
        }
    }
}
