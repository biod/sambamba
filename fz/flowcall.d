module fz.flowcall;

import BioD.Base;

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
}

/// Get flow calls from signal intensities and flow order.
auto flowCalls(ushort[] intensities, string flow_order) {
    
    static FlowCall flowCall(T)(T call) {
        return FlowCall(call[0], Base(call[1]));
    }

    return map!flowCall(zip(intensities, flow_order));
}

/// Get read flow calls. Takes ZF tag and strandness into account.
///
/// Tag name is an optional argument because it is not standard and will likely
/// be changed in the future (there was a proposal on samtools mailing list
/// to introduce standard FB tag).
InputRange!ReadFlowCall readFlowCalls(R)(R read, string flow_order, string tag="ZF") {

    static auto readFlowCall(T)(T tup) {
        auto base = tup[0][0][0];
        auto called_length = tup[0][1];
        auto flow_intensity = tup[1];
        auto offset = tup[2];

        return ReadFlowCall(flow_intensity, cast(ushort)offset, cast(ushort)called_length, base);
    }

    static auto readFlowCallRange(S)(S seq, ushort[] intensities, string flow_order)
    {
        auto flow_index = flowIndex(seq, flow_order);
        auto bases_with_flow_indices = zip(seq, flow_index);
        auto homopolymers = group(bases_with_flow_indices);

        auto flow_intensities = indexed(intensities, uniq(flow_index));
        auto offsets = chain(repeat(0, 1), prefixSum(map!"a[1]"(homopolymers)));

        return map!readFlowCall(zip(homopolymers, flow_intensities, offsets));
    }

    auto zf = cast(int)read[tag];
    auto fz = cast(ushort[])read["FZ"];
    flow_order = flow_order[zf .. $];
    auto intensities = fz[zf .. $];
    if (!read.is_reverse_strand) {
        auto seq = read.sequence;
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order));
    } else {
        auto seq = retro(map!"a.complement"(read.sequence));
        return inputRangeObject(readFlowCallRange(seq, intensities, flow_order));
    }
}
