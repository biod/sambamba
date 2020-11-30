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
module bio.core.genotype;

import bio.core.base;
import bio.core.tinymap;

/// Holds ordered pair of two alleles
struct DiploidGenotype(B) {

    mixin TinyMapInterface!(B.ValueSetSize ^^ 2);

    private static ubyte _getCode(B b1, B b2) {
        auto c1 = b1.internal_code;
        auto c2 = b2.internal_code;
        return cast(ubyte)(c1 * B.ValueSetSize + c2);
    }

    /// Construct a genotype from two bases
    /// Every ambiguous base gets converted to 'N' internally.
    this(B b1, B b2) {
        _code = _getCode(b1, b2);
    }

    /// Construct homozygous genotype
    this(B b) {
        _code = _getCode(b, b);
    }

    /// First allele
    B base1() @property const {
        return B.fromInternalCode(_code / B.ValueSetSize);
    }

    /// Second allele
    B base2() @property const {
        return B.fromInternalCode(_code % B.ValueSetSize);
    }

    ///
    bool is_heterozygous() @property const {
        return base1 != base2;
    }

    ///
    bool is_homozygous() @property const {
        return base1 == base2;
    }

    /// String representation B1|B2 (TODO: add phasing in future?)
    string toString() const {
        return base1 ~ "|" ~ base2;
    }
}

/// Create an instance of DiploidGenotype
auto diploidGenotype(B...)(B bases) {
    return DiploidGenotype!(B[0])(bases);
}

unittest {
    auto g1 = diploidGenotype(Base('C'), Base('W'));
    assert(g1.base1 == 'C');
    assert(g1.base2 == 'W');

    // By default, Base5 is used
    auto g2 = diploidGenotype(Base5('C'));
    assert(g2.base1 == g2.base2);

    // Both bases must be of the same type
    static assert(!__traits(compiles, diploidGenotype(Base5('T'), Base16('D'))));
}
