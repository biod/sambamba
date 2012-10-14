module BioD.Base;

import std.traits;

///
struct Base {
    immutable ValueSetSize = 16;
    static assert(ValueSetSize <= (2 ^^ (_code.sizeof * 8)));

    private ubyte _code;
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

    private static ubyte _getCode(char c) {
        return _char2code[cast(ubyte)c];
    }

    invariant() {
        assert(_code < ValueSetSize);
    }

    /// Construct from IUPAC code
    this(char c) {
        _code = _getCode(c);
    }

    /// Get internal code
    ubyte internal_code() @property const {
        return _code;
    }

    /// Construct from internal code
    static Base fromInternalCode(ubyte code) {
        Base base = void;
        base._code = code;
        return base;
    }

    /// Convert to any kind of char
    T opCast(T)() const 
        if(isSomeChar!T)
    {
        return _code2char[_code];
    }

    /// Convert to char
    char asCharacter() @property const { return opCast!char(); }
    ///
    alias asCharacter this;
}

unittest {
    Base b = 'W';
    assert(b == 'W');

    b = Base.fromInternalCode(0);
    assert(b == '=');
}
