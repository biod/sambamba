/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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

module bio.std.hts.bam.cigar;

import std.algorithm;
import std.range;
import std.conv;
import std.format;
import std.exception;
import std.system;
import std.traits;
import std.array;
import std.bitmanip;
import core.stdc.stdlib;

import bio.core.base;
import bio.core.utils.format;

import bio.std.hts.bam.abstractreader;

import bio.std.hts.bam.writer;
import bio.std.hts.bam.tagvalue;
import bio.std.hts.bam.bai.bin;

import bio.std.hts.bam.md.core;

import bio.std.hts.utils.array;
import bio.std.hts.utils.value;
import bio.core.utils.switchendianness;

import bio.std.hts.thirdparty.msgpack : Packer, unpack;

/**
  Represents single CIGAR operation
 */
struct CigarOperation {
    static assert(CigarOperation.sizeof == uint.sizeof);
    /*
        WARNING!

      It is very essential that the size of
      this struct is EXACTLY equal to uint.sizeof!

      The reason is to avoid copying of arrays during alignment parsing.

      Namely, when some_pointer points to raw cigar data,
      we can just do a cast. This allows to access those data
      directly, not doing any memory allocations.
    */

    private uint raw; // raw data from BAM

    private static ubyte char2op(char c) {
        switch(c) {
            case 'M': return 0;
            case 'I': return 1;
            case 'D': return 2;
            case 'N': return 3;
            case 'S': return 4;
            case 'H': return 5;
            case 'P': return 6;
            case '=': return 7;
            case 'X': return 8;
            default:  return 15; // 15 is used as invalid value
        }
    }

    /// Length must be strictly less than 2^28.
    /// $(BR)
    /// Operation type must be one of M, I, D, N, S, H, P, =, X.
    this(uint length, char operation_type) {
        enforce(length < (1<<28), "Too big length of CIGAR operation");
        raw = (length << 4) | char2op(operation_type);
    }

    this(uint _raw) {
        raw = _raw;
    }

    /// Operation length
    uint length() @property const nothrow @nogc {
        return raw >> 4;
    }

    /// CIGAR operation as one of MIDNSHP=X.
    /// Absent or invalid operation is represented by '?'
    char type() @property const nothrow @nogc {
        return "MIDNSHP=X????????"[raw & 0xF];
    }

    // Each pair of bits has first bit set iff the operation is query consuming,
    // and second bit set iff it is reference consuming.
    //                                            X  =  P  H  S  N  D  I  M
    private static immutable uint CIGAR_TYPE = 0b11_11_00_00_01_10_10_01_11;

    /// True iff operation is one of M, =, X, I, S
    bool is_query_consuming() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 1) != 0;
    }

    /// True iff operation is one of M, =, X, D, N
    bool is_reference_consuming() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 2) != 0;
    }

    /// True iff operation is one of M, =, X
    bool is_match_or_mismatch() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 3) == 3;
    }

    /// True iff operation is one of 'S', 'H'
    bool is_clipping() @property const nothrow @nogc {
        return ((raw & 0xF) >> 1) == 2; // 4 or 5
    }

    void toSam(Sink)(auto ref Sink sink) const
        if (isSomeSink!Sink)
    {
        sink.write(length);
        sink.write(type);
    }

    void toString(scope void delegate(const(char)[]) sink) const {
        toSam(sink);
    }
}

alias CigarOperation[] CigarOperations;

bool is_unavailable(CigarOperations cigar) @property nothrow @nogc {
  return (cigar.length == 1 && cigar[0].raw == '*');
}

/// Forward range of extended CIGAR operations, with =/X instead of M
/// Useful for, e.g., detecting positions of mismatches.
struct ExtendedCigarRange(CigarOpRange, MdOpRange) {
    static assert(isInputRange!CigarOpRange && is(Unqual!(ElementType!CigarOpRange) == CigarOperation));
    static assert(isInputRange!MdOpRange && is(Unqual!(ElementType!MdOpRange) == MdOperation));

    private {
        CigarOpRange _cigar;
        MdOpRange _md_ops;
        CigarOperation _front_cigar_op;
        MdOperation _front_md_op;
        uint _n_mismatches;
        bool _empty;
    }

    ///
    this(CigarOpRange cigar, MdOpRange md_ops) {
        _cigar = cigar;
        _md_ops = md_ops;
        fetchNextCigarOp();
        fetchNextMdOp();
    }

    /// Forward range primitives
    bool empty() @property const {
        return _empty;
    }

    /// ditto
    CigarOperation front() @property {
        debug {
            import std.stdio;
            writeln(_front_cigar_op, " - ", _front_md_op);
        }

        if (_front_cigar_op.type != 'M')
            return _front_cigar_op;

        if (_n_mismatches == 0) {
            assert(_front_md_op.is_match);
            uint len = min(_front_md_op.match, _front_cigar_op.length);
            return CigarOperation(len, '=');
        }

        assert(_front_md_op.is_mismatch);
        return CigarOperation(min(_n_mismatches, _front_cigar_op.length), 'X');
    }

    /// ditto
    ExtendedCigarRange save() @property {
        typeof(return) r = void;
        r._cigar = _cigar.save;
        r._md_ops = _md_ops.save;
        r._front_cigar_op = _front_cigar_op;
        r._front_md_op = _front_md_op;
        r._n_mismatches = _n_mismatches;
        r._empty = _empty;
        return r;
    }

    /// ditto
    void popFront() {
        if (!_front_cigar_op.is_match_or_mismatch) {
            if (_front_cigar_op.is_reference_consuming)
                fetchNextMdOp();
            fetchNextCigarOp();
            return;
        }

        auto len = _front_cigar_op.length;
        if (_n_mismatches > 0) {
            enforce(_front_md_op.is_mismatch);

            if (len > _n_mismatches) {
                _front_cigar_op = CigarOperation(len - _n_mismatches, 'M');
                _n_mismatches = 0;
                fetchNextMdOp();
            } else if (len < _n_mismatches) {
                _n_mismatches -= len;
                fetchNextCigarOp();
            } else {
                fetchNextCigarOp();
                fetchNextMdOp();
            }
        } else {
            enforce(_front_md_op.is_match);
            auto n_matches = _front_md_op.match;

            if (len > n_matches) {
                _front_cigar_op = CigarOperation(len - n_matches, 'M');
                fetchNextMdOp();
            } else if (len < n_matches) {
                _front_md_op.match -= len;
                fetchNextCigarOp();
            } else {
                fetchNextCigarOp();
                fetchNextMdOp();
            }
        }
    }

    private {
        void fetchNextCigarOp() {
            if (_cigar.empty) {
                _empty = true;
                return;
            }

            _front_cigar_op = _cigar.front;
            _cigar.popFront();
        }

        void fetchNextMdOp() {
            if (_md_ops.empty)
                return;

            _n_mismatches = 0;

            _front_md_op = _md_ops.front;
            _md_ops.popFront();

            if (_front_md_op.is_mismatch) {
                _n_mismatches = 1;
                while (!_md_ops.empty && _md_ops.front.is_mismatch) {
                    _md_ops.popFront();
                    _n_mismatches += 1;
                }
            }
        }
    }
}

auto makeExtendedCigar(CigarOpRange, MdOpRange)(CigarOpRange cigar, MdOpRange md_ops) {
    return ExtendedCigarRange!(CigarOpRange, MdOpRange)(cigar, md_ops);
}
