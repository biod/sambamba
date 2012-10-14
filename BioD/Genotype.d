module BioD.Genotype;

import BioD.Base;

///
struct DiploidGenotype {

    immutable ValueSetSize = 25;
    static assert(ValueSetSize <= (2 ^^ (_code.sizeof * 8)));

    private ubyte _code;

    immutable ubyte[16] nt16_to_nt5 = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];
    immutable ubyte[5] nt5_to_nt16 = [1, 2, 4, 8, 15];
    immutable nt5_to_char = "ACGTN";

    private static ubyte _getCode(Base b1, Base b2) {
        auto c1 = nt16_to_nt5[b1.internal_code];
        auto c2 = nt16_to_nt5[b2.internal_code];
        return cast(ubyte)(c1 * nt5_to_char.length + c2);
    }

    /// Construct a genotype from two bases
    /// Every ambiguous base gets converted to 'N' internally.
    this(Base b1, Base b2) {
        _code = _getCode(b1, b2);
    }

    /// Get internal code
    ubyte internal_code() @property const {
        return _code;
    }

    /// Construct from internal code
    static DiploidGenotype fromInternalCode(ubyte code) {
        DiploidGenotype g = void;
        g._code = code;
        return g;
    }

    /// First allele
    Base base1() @property const {
        return Base.fromInternalCode(nt5_to_nt16[_code / nt5_to_char.length]);
    }

    /// Second allele
    Base base2() @property const {
        return Base.fromInternalCode(nt5_to_nt16[_code % nt5_to_char.length]);
    }
}

unittest {
    auto g = DiploidGenotype(Base('C'), Base('W'));
    assert(g.base1 == 'C');
    assert(g.base2 == 'N'); // implicit conversion

    g = DiploidGenotype(Base('A'), Base('A'));
    assert(g.base1 == g.base2);
}
