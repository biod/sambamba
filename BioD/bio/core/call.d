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
module bio.core.call;

import bio.core.base;
import bio.core.genotype;

/// A genotype call
struct Call(alias Gt, B) 
{
    alias Gt!B G;

    private {
        string _sample = void;
        string _chr = void;
        ulong _pos = void;
        B _refbase = void;
        G _gt = void;
        float _qual = void;
    }

    /// Constructor
    this(string sample, string chr, ulong pos,
         B refbase, G genotype, float quality=float.nan)
    {
        _sample = sample;
        _chr = chr;
        _pos = pos;
        _refbase = refbase;
        _gt = genotype;
        _qual = quality;
    }

    /// Sample name
    string sample() @property const {
        return _sample;
    }

    /// Chromosome name
    string chromosome() @property const {
        return _chr;
    }

    /// 0-based position on the reference
    ulong position() @property const {
        return _pos;
    }

    /// Reference base at the site
    B reference_base() @property const {
        return _refbase;
    }

    /// Most probable genotype
    ref const(G) genotype() @property const {
        return _gt;
    }

    ///
    alias genotype this; 

    /// Phred-scaled quality score. If unknown, set to NaN.
    float quality() @property const {
        return _qual;
    }

    /// Returns true if this call is not a reference one.
    bool is_variant() @property const {
        return _gt != G(_refbase);
    }
}

alias Call!(DiploidGenotype, Base5) DiploidCall5;
alias Call!(DiploidGenotype, Base16) DiploidCall;
alias DiploidCall DiploidCall16;

unittest {
    auto call = DiploidCall("NA01234", "chr10", 543210,
                            Base('T'), diploidGenotype(Base('C'), Base('T')),
                            47.0);

    assert(call.is_variant);
    assert(call.is_heterozygous);
    assert(call.reference_base == 'T');
}
