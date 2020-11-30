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
module bio.core.base;

import bio.core.tinymap;
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

    private enum ubyte[256] _char2code = [
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
    private enum _code2char = "=ACMGRSVTWYHKDBN";

    private enum ubyte[16] _complement_table = [0x0, 0x8, 0x4, 0xC, 
                                                     0x2, 0xA, 0x6, 0xE,
                                                     0x1, 0x9, 0x5, 0xD,
                                                     0x3, 0xB, 0x7, 0xF];
    /// Complementary base
    Base complement() @property const {
        // take the code, reverse the bits, and return the base
        return Base.fromInternalCode(_complement_table[_code]);
    }

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

    /// ditto
    this(dchar c) {
        _code = _char2code[cast(ubyte)c];
    }

    private enum ubyte[5] nt5_to_nt16 = [1, 2, 4, 8, 15];
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

    T opCast(T)() const 
        if (is(Unqual!T == char) || is(Unqual!T == dchar))
    {
        return asCharacter;
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
/// (internal codes are 0, 1, 2, 3, and 4 correspondingly)
struct Base5 {
    mixin TinyMapInterface!5;

    private enum ubyte[256] _char2code = [
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

    private enum _code2char = "ACGTN";
    private enum ubyte[16] nt16_to_nt5 = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

    mixin CommonBaseOperations;

    /// Complementary base
    Base5 complement() @property const {
        return Base5.fromInternalCode(cast(ubyte)(_code == 4 ? 4 : (3 - _code)));
    }

    /// Construct base from one of "acgtACGT" symbols.
    /// Every other character is converted to 'N'
    this(char c) {
        _code = _char2code[cast(ubyte)c];
    }

    /// ditto
    this(dchar c) {
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

    T opCast(T)() const 
        if (is(Unqual!T == char) || is(Unqual!T == dchar))
    {
        return asCharacter;
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

/// Complement base, which might be Base5, Base16, char, or dchar.
B complementBase(B)(B base) {
    static if(is(Unqual!B == dchar) || is(Unqual!B == char))
    {
        return cast(B)(Base16(base).complement);
    }
    else
        return base.complement;
}

/// Convert character to base
template charToBase(B=Base16)
{
    B charToBase(C)(C c)
        if(is(Unqual!C == char) || is(Unqual!C == dchar))
    {
        return B(c);
    }
}

unittest {
    assert(complementBase('T') == 'A');
    assert(complementBase('G') == 'C');

    assert(complementBase(Base5('A')) == Base5('T'));
    assert(complementBase(Base16('C')) == Base16('G'));

    assert(charToBase!Base16('A').complement == Base16('T'));
}
