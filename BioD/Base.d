module BioD.Base;

import BioD.TinyMap;
import std.traits;

/// Code common to both Base5 and Base16
mixin template CommonBaseOperations() {
    /// Convert to char
    char asCharacter() @property const { return _code2char[_code]; }
    ///
    alias asCharacter this;
}

/// Base representation supporting full set of IUPAC codes
struct Base {
    mixin TinyMapInterface!16; 

    private immutable ubyte[256] _char2code = [
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,

                15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
                15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
                15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
                15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,

                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
            ];

    // = 0000
    //
    // A 0001
    // C 0010
    // G 0100
    // T 1000
    //
    // W 1001 (A T) Weak
    // S 0110 (C G) Strong
    //
    // M 0011 (A C) aMino
    // K 1100 (G T) Keto
    // R 0101 (A G) puRine
    // Y 1010 (A G) pYrimidine
    //
    // B 1110 (not A)
    // D 1101 (not C)
    // H 1011 (not G)
    // V 0111 (not T)
    //
    // N 1111 (aNy base)
    private immutable _code2char = "=ACMGRSVTWYHKDBN";

    unittest {
        import std.ascii;

        foreach (i, c; _code2char) {
            assert(_char2code[c] == i);
        }

        foreach (c; 0 .. 256) {
            auto c2 = _code2char[_char2code[c]];
            if (c2 != 'N') {
                if ('0' <= c && c <= '9') {
                    assert(c2 == "ACGT"[c - '0']);
                } else {
                    assert(c2 == toUpper(c));
                }
            }
        }
    }

    mixin CommonBaseOperations;
    /// Construct from IUPAC code
    this(char c) {
        _code = _char2code[cast(ubyte)c];
    }

    private immutable ubyte[5] nt5_to_nt16 = [1, 2, 4, 8, 15];
    private static Base fromBase5(Base5 base) {
        Base b = void;
        b._code = nt5_to_nt16[base.internal_code];
        return b;
    }

    /// Conversion to Base5
    Base5 opCast(T)() const
        if (is(T == Base5)) 
    {
        return Base5.fromBase16(this);
    }
}

unittest {
    Base b = 'W';
    assert(b == 'W');

    b = Base.fromInternalCode(0);
    assert(b == '=');
}

alias Base Base16;

/// Base representation supporting only 'A', 'C', 'G', 'T', and 'N'
struct Base5 {
    mixin TinyMapInterface!5;

    private immutable ubyte[256] _char2code = [
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
                4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
                ];

    private immutable _code2char = "ACGTN";
    private immutable ubyte[16] nt16_to_nt5 = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

    mixin CommonBaseOperations;

    /// Construct base from one of "acgtACGT" symbols.
    /// Every other character is converted to 'N'
    this(char c) {
        _code = _char2code[cast(ubyte)c];
    }

    private static Base5 fromBase16(Base16 base) {
        Base5 b = void;
        b._code = nt16_to_nt5[base.internal_code];
        return b;
    }

    /// Conversion to Base16
    Base16 opCast(T)() const
        if(is(T == Base16)) 
    {
        return Base16.fromBase5(this);
    }
}

unittest {
    auto b5 = Base5('C');
    assert(b5.internal_code == 1);
    b5 = Base5.fromInternalCode(3);
    assert(b5 == 'T');

    // doesn't work with std.conv.to
    //
    //import std.conv;
    //assert(to!Base16(b5).internal_code == 8);

    assert((cast(Base16)b5).internal_code == 8);
}
