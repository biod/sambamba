/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

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
module bio.std.hts.iontorrent.flowcall;
import bio.std.hts.bam.tagvalue;
import bio.std.hts.iontorrent.flowindex;

import bio.core.base;
import bio.core.utils.range;

import std.array;
import std.typecons;
import std.range;
import std.algorithm;
import std.exception;

/// Tag where flow signal intensities are stored
enum FlowGramTag : ubyte {
    FZ,
    ZM
}

/// Scale of intensity values
float multiplier(FlowGramTag tag) {
    return tag == FlowGramTag.FZ ? 100.0 : 256.0;
}

/// Flow base call
struct FlowCall {
    private {
        short _signal_intensity;

        static assert(Base.ValueSetSize <= 16 && FlowGramTag.max < 16,
                      "implementation of FlowCall should be changed?");

        ubyte _storage;  // tag in upper 4 bits, base in lower 4 bits

        Base _base() @property const {
            return Base.fromInternalCode(_storage & 0xF);
        }

        void _base(Base b) @property {
            _storage &= 0xF0;
            _storage |= b.internal_code;
        }

        FlowGramTag _tag() @property const {
            return cast(FlowGramTag)(_storage >> 4);
        }

        void _tag(FlowGramTag tag) @property {
            _storage &= 0xF;
            _storage |= (cast(ubyte)tag << 4);
        }

        this(short signal_intensity, Base b, FlowGramTag tag) {
            _signal_intensity = signal_intensity;
            _storage = cast(ubyte)(b.internal_code | (tag << 4));
        }
    }

    /// Nucleotide
    Base base() @property const {
        return _base;
    }

    /// Signal intensity, normalized to homopolymer lengths
    float intensity() @property const {
        return _signal_intensity / multiplier(_tag);
    }

    /// round(intensity * Multiplier) where Multiplier is 100.0 for FZ tag,
    /// and 256.0 for ZM tag.
    /// More efficient, because this is how intensities are stored in FZ/ZM tag.
    short intensity_value() @property const {
        return _signal_intensity;
    }
}

/// Flow call associated with a read
struct ReadFlowCall {
    private {
        FlowCall _fc;
        ushort _offset;
        ushort _called_len;
        ushort _flow_index;

        this(Base b, short signal_intensity, ushort offset,
                     ushort called, ushort flow_index, FlowGramTag tag)
        {
            _fc = FlowCall(signal_intensity, b, tag);
            _offset = offset;
            _called_len = called;
            _flow_index = flow_index;
        }
    }

    /// Called nucleotide
    Base base() @property const {
        return _fc._base;
    }

    /// Set base to its complement
    void complement() {
        _fc._base = _fc._base.complement;
    }

    /// Called homopolymer length
    ushort length() @property const {
        return _called_len;
    }

    /// Zero-based position of the first nucleotide in the run,      
    /// relative to start of the read. Takes strandness into account.
    ushort offset() @property const {
        return _offset;
    }

    /// Signal intensity, normalized to homopolymer lengths
    float intensity() @property const {
        return _fc.intensity;
    }

    /// round(intensity * Multiplier) where Multiplier is 100.0 for FZ tags,
    /// and 256.0 for ZM tags.
    /// More efficient, because this is how intensities are stored in FZ/ZM tag.
    short intensity_value() @property const {
        return _fc._signal_intensity;
    }

    /// Flow index (0-based)
    size_t flow_index() @property const {
        return _flow_index;
    }
}

/// Get flow calls from signal intensities and flow order.
auto flowCalls(short[] intensities, string flow_order, FlowGramTag tag) {
    
    static FlowCall flowCall(T)(T call) {
        return FlowCall(call[0], Base(call[1]), call[2]);
    }

    return map!flowCall(zip(intensities, flow_order, repeat(tag)));
}

struct ReadFlowCallRange(S) 
    if (!is(S == class))
{
    private {
        string _flow_order = void;
        short[] _intensities = void;
        bool _rev = void;
        S _sequence = void;

        int _zf = void;
        Base _current_base = void;
        ushort _current_length = void;
        size_t _current_flow_index;
        ushort _current_offset;

        ushort _overlap = void;
        
        FlowGramTag _tag = void;

        bool _empty = false;

        // consumes next homopolymer from the sequence,
        // and updates _current_base, _current_flow_index, 
        // _current_length appropriately
        void _doSetup() {
            if (_sequence.empty) {
                _empty = true;
                return;
            }

            _current_length = 1; 

            // setup current base and current length
            if (!_rev) {
                _current_base = _sequence.front;
                _sequence.popFront();
                while (!_sequence.empty && _sequence.front == _current_base) {
                    _sequence.popFront();
                    ++_current_length;
                }
            } else {
                _current_base = _sequence.back; // complement later
                _sequence.popBack();            // because of comparison below
                while (!_sequence.empty && _sequence.back == _current_base) {
                    _sequence.popBack();
                    ++_current_length;
                }
                _current_base = _current_base.complement;
            }

            // setup current flow index
            for ( ; _current_flow_index < _flow_order.length; ++_current_flow_index) {
                if (_flow_order[_current_flow_index] == _current_base) {
                    break;
                }
            }
        }
    }

    this(S seq, short[] intensities, bool reverse_strand, 
         string flow_order, ushort first_base_overlap, int zf, FlowGramTag tag) 
    {
        _sequence = seq;
        _intensities = intensities;
        _rev = reverse_strand;
        _flow_order = flow_order;
        _zf = zf;
        _overlap = first_base_overlap;
        _tag = tag;

        if (_sequence.empty) {
            _empty = true;
        } else {
            _doSetup();
        }
    }

    bool empty() @property const {
        return _empty;
    }

    ReadFlowCall front() @property const {
        enforce(_current_flow_index < _intensities.length,
                "Inconsistency between FZ/ZM tag and read bases");

        auto intensity = cast(ushort)(_intensities[_current_flow_index] - _overlap);
        ReadFlowCall rfc = void;
        rfc._fc = FlowCall(intensity, _current_base, _tag);
        rfc._offset = _current_offset;
        rfc._called_len = _current_length;
        rfc._flow_index = cast(ushort)(_current_flow_index + _zf);
        return rfc;
    }

    void popFront() {
        _current_offset += _current_length;

        ++_current_flow_index;
        _overlap = 0; // after first base it is always zero

        _doSetup();
    }

    ReadFlowCallRange!S save() @property {
        // bitwise copy
        // FIXME: is it safe?
        ReadFlowCallRange!S r = this;
        return r;
    }
}

private ReadFlowCallRange!S readFlowCallRange(S)(S seq, short[] intensities, bool rev,
                                                 string flow_order, ushort overlap, int zf,
                                                 FlowGramTag tag)
{
    return ReadFlowCallRange!S(seq, intensities, rev, flow_order, overlap, zf, tag);
}


/// Get read flow calls. Takes ZF tag and strandness into account.
///
/// Tag name is an optional argument because it is not standard and will likely
/// be changed in the future (there was a proposal on samtools mailing list
/// to introduce standard FB tag).
auto readFlowCalls(R)(R read, string flow_order, string key_sequence, string tag="ZF") {

    auto zf = cast(int)read[tag];
    auto fz_value = read["FZ"];
    auto zm_value = read["ZM"];

    enforce(!(fz_value.is_nothing && zm_value.is_nothing),
            "Neither FZ nor ZM tag is presented in a mapped read");

    auto fg_tag = fz_value.is_nothing ? FlowGramTag.ZM : FlowGramTag.FZ;

    short[] flow_int = *cast(short[]*)(fg_tag == FlowGramTag.ZM ? &zm_value : &fz_value);

    flow_order = flow_order[zf .. $];
    auto intensities = flow_int[zf .. $];

    // key sequence is required because its last base can overlap with first called base
    ushort overlap = 0;

    Base5 base = read.is_reverse_strand ? read.sequence.back.complement : read.sequence.front;
    foreach_reverse (c; key_sequence) {
        if (c != base)
            break;
        overlap += cast(int)(multiplier(fg_tag));
    }

    return readFlowCallRange(read.sequence, intensities, read.is_reverse_strand,
                             flow_order, overlap, zf, fg_tag);
}
