/*
    This file is part of BioD.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.sam.reader;
//
import bio.std.hts.bam.abstractreader;
import bio.std.hts.sam.header;
import bio.std.hts.bam.read;
import bio.std.hts.bam.reference;
import bio.std.hts.bam.referenceinfo;
import bio.core.utils.outbuffer;
import bio.core.utils.range;

import bio.core.utils.bylinefast;
alias ByLineFast _LineRange;

version(DigitalMars) {
    import bio.std.hts.sam.utils.recordparser;
} else {
    import bio.std.hts.sam.utils.fastrecordparser;
}

import std.stdio;
import std.array;
import std.string;
import std.range;
import std.algorithm;
import std.typecons;
import std.parallelism;
import std.process;
import std.exception;
import core.stdc.string;

BamRead _parseSamRecord(Tuple!(char[], SamReader, OutBuffer) t) {
    auto r = parseAlignmentLine(cast(string)t[0], t[1]._header, t[2]);
    BamRead result;
    if (t[1]._seqprocmode) {
        result = r;
    } else {
        auto storage = uninitializedArray!(ubyte[])(r.raw_data.length);
        storage[] = r.raw_data[];
        result.raw_data = storage;
    }
    result.associateWithReader(t[1]);
    return result;
}

private {
    extern(C) size_t lseek(int, size_t, int);
    bool isSeekable(ref File file) {
        return lseek(file.fileno(), 0, 0) != ~0;
    }
}

///
class SamReader : IBamSamReader {

    private {
        version(gzippedSamSupport) {
        void checkGunzip() {
            auto gunzip = executeShell("gunzip -V");
            if (gunzip.status != 0)
                throw new Exception("gunzip is not installed on this system, can't read gzipped SAM");
        }

        File openSamFile(string filename) {
            if (filename.length < 4)
                throw new Exception("invalid name for SAM file: " ~ filename);
            if (filename[$ - 3 .. $] == ".gz") {
                checkGunzip();
                auto pipe = pipeShell("gunzip -c " ~ filename);
                return pipe.stdout;
            } else if (filename[$ - 4 .. $] == ".bam") {
                throw new Exception("SAM reader can't read BAM file " ~ filename);
            } else {
                return File(filename);
            }
        }

        } else {

        File openSamFile(string filename) {
            if (filename[$ - 4 .. $] == ".bam") {
                throw new Exception("SAM reader can't read BAM file " ~ filename);
            } else {
                return File(filename);
            }
        }

        }
    }

    ///
    this(string filename) {
        _file = openSamFile(filename);
        _filename = filename;
        _seekable = _file.isSeekable();
        _initializeStream();
    }

    ///
    bio.std.hts.sam.header.SamHeader header() @property {
        return _header;
    }

    ///
    const(bio.std.hts.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences() @property const {
        return _reference_sequences;
    }

    ///
    bool hasReference(string reference) {
        return null != (reference in _reference_sequence_dict);
    }

    ///
    bio.std.hts.bam.reference.ReferenceSequence opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " is not present in the header");
        auto ref_id = _reference_sequence_dict[ref_name];
        return ReferenceSequence(null, ref_id, _reference_sequences[ref_id]);
    }

    /// Reads in SAM file.
    auto reads() @property {

        _LineRange lines = _lines;
        if (_seekable) {
            if (_filename !is null) {
                auto file = openSamFile(_filename);
                lines = ByLineFast(file);
            } else {
                _file.seek(0);
                lines = ByLineFast(_file);
            }
            auto dummy = lines.front;
            for (int i = 0; i < _lines_to_skip; i++)
                lines.popFront();
        }

        auto b = new OutBuffer(262144);
        return lines.zip(repeat(this), repeat(b)).map!_parseSamRecord();
    }

    ///
    void assumeSequentialProcessing() {
        _seqprocmode = true;
    }

    ///
    std.range.InputRange!(bio.std.hts.bam.read.BamRead) allReads() @property {
        return inputRangeObject(reads);
    }

    /// Filename
    string filename() @property const {
        return _filename;
    }
private:

    File _file;
    bool _seekable;
    string _filename;
    _LineRange _lines;
    ulong _lines_to_skip;

    bool _seqprocmode;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict;

    void _initializeStream() {
        auto header = Appender!(char[])();

        _lines = ByLineFast(_file);

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
                header.put('\n');
                _lines_to_skip += 1;
                _lines.popFront();
            } else {
                break;
            }
        }

        import core.memory;
        GC.disable();
        _header = new SamHeader(cast(string)(header.data));
        GC.enable();

        _reference_sequences = new ReferenceSequenceInfo[_header.sequences.length];
        foreach (sq; _header.sequences) {
            auto seq = ReferenceSequenceInfo(sq.name, sq.length);
            auto n = cast(int)_reference_sequences.length;
            _reference_sequence_dict[sq.name] = n;
            _reference_sequences[_header.getSequenceIndex(seq.name)] = seq;
        }
    }
}
