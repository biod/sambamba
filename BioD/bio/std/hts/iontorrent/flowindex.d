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
module bio.std.hts.iontorrent.flowindex;

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
    import bio.core.base;

    import std.conv;
    import std.algorithm;

    auto seq = map!(c => Base5(c))("AACGTAAACCTCACT");
    string flow_order = "ATGCATGCATGCATGCATGCATGCATGC";
                      // 0123456789111111111122222222
                      //           012345678901234567
    assert(equal(flowIndex(seq, flow_order), [0, 0, 3, 6, 9, 12, 12, 12, 15, 15, 17, 19, 20, 23, 25]));
}
