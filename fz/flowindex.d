module fz.flowindex;

import std.range;

///
struct FlowIndex(S) 
{
    private {
        S _seq;
        string _fo;
        size_t _index;
    }

    this(S sequence, string flow_order) {
        _seq = sequence;
        _fo = flow_order;

        if (!_seq.empty) {
            while (_index < _fo.length) {
                if (_fo[_index] == _seq.front)
                    break;

                ++_index;
            }
        }
    }

    ///
    bool empty() @property {
        return _seq.empty || (_index == _fo.length);
    }

    /// Current flow index
    size_t front() @property const {
        return _index;
    }

    /// Move to next read base
    void popFront() {
        auto prev_base = _seq.front;
        _seq.popFront();
        if (_seq.empty) return;

        if (_seq.front == prev_base) {
            return; // keep index as is
        }

        _index += 1;
        while (_index < _fo.length) {
            if (_fo[_index] == _seq.front)
                break;

            _index++;
        }
    }
}

/// Given a sequence of bases and flow order, recover flow index,
/// i.e. sequence of 0-based flow positions for each base.
auto flowIndex(S)(S sequence, string flow_order) 
{
    return FlowIndex!S(sequence, flow_order);
}

unittest {
    import BioD.Base;

    import std.conv;
    import std.algorithm;

    auto seq = map!(c => Base5(c))("AACGTAAACCTCACT");
    string flow_order = "ATGCATGCATGCATGCATGCATGCATGC";
                      // 0123456789111111111122222222
                      //           012345678901234567
    assert(equal(flowIndex(seq, flow_order), [0, 0, 3, 6, 9, 12, 12, 12, 15, 15, 17, 19, 20, 23, 25]));
}
