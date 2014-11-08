module cram.reader;

import cram.htslib;
import cram.readrange;
import cram.reference;
import cram.exception;

import bio.bam.abstractreader, bio.bam.referenceinfo,
       bio.sam.header, bio.bam.reference, bio.bam.read;
import std.string, std.array, std.range, std.conv;

auto zeroChecked(alias func, T...)(string err_msg, auto ref T params) {
    int ret = func(params);
    if (!ret)
        throw new CramException(err_msg);
}

class CramReader : IBamSamReader {
    private {
        string _fn;
        string _mode;
        cram_fd* _fd;
        int _n_threads;

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

    private auto openCram(string filename) {
        auto fd = cram_open(toStringz(filename), "rb");
        cram_set_option(fd, cram_option.CRAM_OPT_DECODE_MD);
// the library is not ready yet for this :C
//        if (_n_threads > 1)
//          cram_set_option(fd, cram_option.CRAM_OPT_NTHREADS, _n_threads);
        return fd;
    }

    this(string filename, uint n_threads) {
        _fn = filename;
        _n_threads = n_threads;

        _fd = openCram(filename);

        if (_fd == null) {
            throw new CramException("can't open file " ~ filename);
        }

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
        return cram.reference.ReferenceSequence(this, fd(), _seq_op, ref_id,
                                                _reference_sequences[ref_id]);
    }

    ///
    cram.reference.ReferenceSequence reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return cram.reference.ReferenceSequence(this, fd(), _seq_op, ref_id,
                                                _reference_sequences[ref_id]);
    }

    BamRead[] unmappedReads() {
        throw new Exception("* region unimplemented for CRAM files");
        return [];
    }

    SamHeader header() @property {
        if (_header is null)
            _header = parseHeader();
        return _header;
    }

    const(ReferenceSequenceInfo)[] reference_sequences() @property const {
        return _reference_sequences;
    }

    private bool _seq_op;
    void assumeSequentialProcessing() {
        _seq_op = true;
    }

    void close() {
        zeroChecked!cram_close("failed to close CRAM file descriptor",
                _fd);
    }

    private cram_fd* fd() {
        return _fn == "-" ? _fd : openCram(_fn);
    }

    CramBamReadRange reads() {
        return typeof(return)(fd(), _seq_op, this);
    }

    std.range.InputRange!(bio.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads());
    }

    void createIndex() {
        int ret = cram_index_build(_fd, toStringz(_fn));
        if (ret != 0) {
            throw new Exception("failed to build index for CRAM file " ~ _fn);
        }
    }
}
