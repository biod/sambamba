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
module bio.core.tinymap;

private import std.algorithm;
private import std.range;
private import std.traits;

import std.bitmanip;

/// Efficient dictionary for cases when the number of possible keys is small 
/// and is known at compile-time. The data is held in a static array, no allocations occur.
///
/// Key type must:
///     have static member ValueSetSize (integer)
///     have static method fromInternalCode (returning an instance of Key type);
///     have property $(D internal_code) that maps it to integers 0, 1, ..., ValueSetSize - 1
struct TinyMap(K, V, alias TinyMapPolicy=useBitArray) {
    private V[K.ValueSetSize] _dict;
    private size_t _size;
    private mixin TinyMapPolicy!(K, V) Policy;
    private alias ReturnType!(K.internal_code) TCode;

    /// Constructor
    static TinyMap!(K, V, TinyMapPolicy) opCall(Args...)(Args args) {
        TinyMap!(K, V, TinyMapPolicy) result;
        result.Policy.init(args);
        return result;
    }
   
    /// Current number of elements
    size_t length() @property const {
        return _size;
    }

    /// Indexed access
    auto ref V opIndex(Key)(auto ref Key key)
        if(is(Unqual!Key == K))
    {
        assert(key in this);
        return _dict[key.internal_code];
    }

    /// ditto
    auto ref const(V) opIndex(Key)(auto ref Key key) const
        if(is(Unqual!Key == K))
    {
        assert(key in this);
        return _dict[key.internal_code];
    }


    /// ditto
    V opIndexAssign(V value, K key) {
        if (key !in this) {
            ++_size;
        }
        _dict[key.internal_code] = value;
        Policy._onInsert(key);
        return value;
    }

    /// ditto
    void opIndexOpAssign(string op)(V value, K key) {
        if (key !in this) {
            ++_size;
            _dict[key.internal_code] = V.init;
        }
        mixin("_dict[key.internal_code] " ~ op ~ "= value;");
        Policy._onInsert(key);
    }

    /// Check if the key is in the dictionary
    bool opIn_r(K key) const {
        return Policy._hasKey(key);
    }

    /// Removal
    bool remove(K key) {
        if (key in this) {
            --_size;
            Policy._onRemove(key);
            return true;
        }
        return false;
    }

    /// Range of keys
    auto keys() @property const {
        // FIXME: create nice workaround for LDC bug #217
        K[] _ks;
        foreach (i; 0 .. K.ValueSetSize) {
            if (Policy._hasKeyWithCode(i))
                _ks ~= K.fromInternalCode(cast(TCode)i);
        }
        return _ks;
    }

    /// Range of values
    auto values() @property const {
        V[] _vs;
        foreach (i; 0 .. K.ValueSetSize) {
            if (Policy._hasKeyWithCode(i))
                _vs ~= _dict[i];
        }
        return _vs;
    }

    /// Iteration with foreach
    int opApply(scope int delegate(V value) dg) {
        foreach (i; iota(K.ValueSetSize)) {
            if (Policy._hasKeyWithCode(i)) {
                auto ret = dg(_dict[i]);
                if (ret != 0) return ret;
            }
        }
        return 0;
    }

    /// ditto
    int opApply(scope int delegate(K key, V value) dg) {
        foreach (i; iota(K.ValueSetSize)) {
            if (Policy._hasKeyWithCode(i)) {
                auto ret = dg(K.fromInternalCode(cast(TCode)i), _dict[i]);
                if (ret != 0) return ret;
            }
        }
        return 0;
    }
}

/// For each possible key store 0 if it's absent in the dictionary,
/// or 1 otherwise. Bit array is used for compactness.
///
/// This is the default option. In this case, size of dictionary is
/// roughly (V.sizeof + 1/8) * K.ValueSetSize
mixin template useBitArray(K, V) {
    private BitArray _value_is_set;

    private void init() {
        _value_is_set.length = K.ValueSetSize;
    }

    private bool _hasKey(K key) const {
        return _value_is_set[key.internal_code];
    }

    private bool _hasKeyWithCode(size_t code) const {
        return _value_is_set[code];
    }

    private void _onInsert(K key) {
        _value_is_set[key.internal_code] = true;
    }

    private void _onRemove(K key) {
        _value_is_set[key.internal_code] = false;
    }
}

/// Use default value specified at construction as an indicator
/// of key absence.
/// That allows to save K.ValueSetSize bits of memory.
///
/// E.g., you might want to use -1 as such indicator if non-negative
/// numbers are stored in the dictionary.
mixin template useDefaultValue(K, V) {
    private V _default_value;

    private void init(V value) {
        _default_value = value;
        if (_default_value != V.init) {
            _dict[] = _default_value;
        }
    }

    private bool _hasKey(K key) const {
        return _dict[key.internal_code] != _default_value;
    }

    private bool _hasKeyWithCode(size_t code) const {
        return _dict[code] != _default_value;
    }

    private void _onInsert(K key) {}

    private void _onRemove(K key) {
        this[key] = _default_value;
    }
}

/// Allows to set up a dictionary which is always full.
mixin template fillNoRemove(K, V) {

    private void init() {
        _size = K.ValueSetSize;
    }

    private void init(V value) {
        _size = K.ValueSetSize;

        for (size_t i = 0; i < _size; ++i)
            _dict[i] = value;
    }

    private bool _hasKey(K key) const {
        return true;
    }

    private bool _hasKeyWithCode(size_t code) const {
        return true;
    }

    private void _onInsert(K key) {}

    private void _onRemove(K key) {
        ++_size;
    }
}

unittest {

    import std.array;
    import bio.core.base;

    void test(M)(ref M dict) {
        auto b1 = Base('A');
        auto b2 = Base('C');
        auto b3 = Base('G');
        auto b4 = Base('T');
        dict[b1] = 2;
        dict[b2] = 3;
        assert(dict.length == 2);
        assert(dict[b1] == 2);
        assert(b2 in dict);
        assert(b3 !in dict);
        assert(b4 !in dict);
        dict[b4] = 5;
        assert(equal(sort(array(dict.values)), [2, 3, 5]));
        dict.remove(b1);
        assert(b1 !in dict);
        assert(dict.length == 2);
        assert(dict[b2] == 3);

        foreach (k, v; dict) {
            assert(k in dict);
            assert(dict[k] == v);
        }
    }

    auto dict1 = TinyMap!(Base, int)();
    auto dict2 = TinyMap!(Base, int, useDefaultValue)(-1);
    int[Base] dict3;

    test(dict1);
    test(dict2);
    test(dict3);

    auto dict4 = TinyMap!(Base, ulong[4])();
    dict4[Base('A')] = [0, 1, 2, 3];
    dict4[Base('A')][3] += 1;
    assert(dict4[Base('A')] == [0, 1, 2, 4]);
}

/// Convenient mixin template for getting your struct working with TinyMap.
///
/// Creates
///     1) private member of type T with name _code
///     2) fromInternalCode static method
///     3) internal_code property
///     4) static member ValueSetSize equal to N
///     5) invariant that _code is always less than ValueSetSize
///
/// That is, the only thing which implementation is up to you is
/// setting _code appropriately.
mixin template TinyMapInterface(uint N, T=ubyte) if (isUnsigned!T) {
    private T _code;

    enum ValueSetSize = N;
    static assert(N <= 2 ^^ (T.sizeof * 8));

    static typeof(this) fromInternalCode(T code) {
        typeof(this) obj = void;
        obj._code = code;
        return obj;
    }

    T internal_code() @property const {
        return _code;
    }
}
