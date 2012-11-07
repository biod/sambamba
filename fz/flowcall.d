module fz.flowcall;

import BioD.Base;

import tagvalue;
import fz.flowindex;
import utils.range;

import std.array;
import std.typecons;
import std.range;
import std.algorithm;

/// Flow base call
struct FlowCall {
    private {
        ushort _signal_intensity;
        Base _base;
    }

    /// Nucleotide
    Base base() @property const {
        return _base;
    }

    /// Signal intensity, normalized to homopolymer lengths
    float intensity() @property const {
        return _signal_intensity / 100.0;
    }

    /// round(intensity * 100.0)
    /// More efficient, because this is how intensities are stored in FZ tag.
    ushort intensity_value() @property const {
        return _signal_intensity;
    }
}

/// Flow call associated with a read
struct ReadFlowCall {
    private {
        ushort _signal_intensity;
        ushort _offset;
        ushort _called_len;
        Base _base;
        size_t _flow_index;
    }

    /// Called nucleotide
    Base base() @property const {
        return _base;
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
        return _signal_intensity / 100.0;
    }

    /// round(intensity * 100.0)
    /// More efficient, because this is how intensities are stored in FZ tag.
    ushort intensity_value() @property const {
        return _signal_intensity;
    }

    /// Flow index (0-based)
    size_t flow_index() @property const {
        return _flow_index;
    }
}

/// Get flow calls from signal intensities and flow order.
auto flowCalls(ushort[] intensities, string flow_order) {
    
    static FlowCall flowCall(T)(T call) {
        return FlowCall(call[0], Base(call[1]));
    }

    return map!flowCall(zip(intensities, flow_order));
}

struct ReadFlowCallRange(S) 
    if (!is(S == class))
{
    private {
        string _flow_order = void;
        ushort[] _intensities = void;
        S _sequence = void;

        int _zf = void;
        Base _current_base = void;
        ushort _current_length = void;
        size_t _current_flow_index;
        ushort _current_offset;

        bool _empty;

        // consumes next homopolymer from the sequence,
        // and updates _current_base, _current_flow_index, 
        // _current_length appropriately
        void _doSetup() {
            if (_sequence.empty) {
                _empty = true;
                return;
            }

            // setup current base and current length
            _current_base = _sequence.front;
            _sequence.popFront();
            _current_length = 1;
            while (!_sequence.empty && _sequence.front == _current_base) {
                _sequence.popFront();
                ++_current_length;
            }

            // setup current flow index
            for ( ; _current_flow_index < _flow_order.length; ++_current_flow_index) {
                if (_flow_order[_current_flow_index] == _current_base) {
                    break;
                }
            }
        }
    }

    this(S seq, ushort[] intensities, string flow_order, int zf) {
        _sequence = seq;
        _intensities = intensities;
        _flow_order = flow_order;
        _zf = zf;

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
        return ReadFlowCall(_intensities[_current_flow_index], _current_offset,
                            _current_length, _current_base, _current_flow_index + _zf);
    }

    void popFront() {
        _current_offset += _current_length;

        ++_current_flow_index;

        _doSetup();
    }

    ReadFlowCallRange!S save() @property {
        // bitwise copy
        // FIXME: is it safe?
        ReadFlowCallRange!S r = this;
        return r;
    }
}

ReadFlowCallRange!S readFlowCallRange(S)(S seq, ushort[] intensities, string flow_order, int zf)
{
    return ReadFlowCallRange!S(seq, intensities, flow_order, zf);
}


/// Get read flow calls. Takes ZF tag and strandness into account.
///
/// Tag name is an optional argument because it is not standard and will likely
/// be changed in the future (there was a proposal on samtools mailing list
/// to introduce standard FB tag).
ForwardRange!ReadFlowCall readFlowCalls(R)(R read, string flow_order, string tag="ZF") {

    static auto readFlowCall(T)(T tup) {
        auto base = tup[0][0][0];
        auto called_length = tup[0][1];
        auto flow_intensity = tup[1];
        auto offset = tup[2];

        return ReadFlowCall(flow_intensity, cast(ushort)offset, cast(ushort)called_length, base);
    }

    auto zf = cast(int)read[tag];
    Value fz_value = read["FZ"];
    ushort[] fz = *cast(ushort[]*)(&fz_value);

    flow_order = flow_order[zf .. $];
    auto intensities = fz[zf .. $];
    if (!read.is_reverse_strand) {
        auto seq = read.sequence;
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order, zf));
    } else {
        auto seq = retro(map!"a.complement"(read.sequence));
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order, zf));
    }
}
