/**
   Contains classes implementing TagStorage interface for accessing tags. 
   The classes encapsulate tag parsing logic as well.

   The approach to parsing used by original samtools is as follows:

   1) don't parse tags eagerly
   2) store ubyte[] array in the struct
   3) change endianness in-place, if needed
        - thus we can slice strings and arrays -- no extra allocations!
   4) to retrieve value, iterate over all tags
        - micro-optimization: compare ushorts instead of 2-char strings


   TODO: consider getting rid of interface for the sake of optimization,
         and make LazyTagStorage a struct.
*/
module tagstorage;

import tagvalue;
import utils.switchendianness;

import std.stream;
import std.conv;
import std.system;

/**
  Abstract tag storage. 

  Provides hash-like access (currently read-only) and opportunity to iterate
  storage like an associative array.
*/
abstract class TagStorage {
    abstract Value opIndex(string s); /// hash-like access

    abstract int opApply(int delegate(ref string key, ref Value) dg); /// iteration
}

import std.system;
import core.exception;

/**
    Provides behaviour similar to original samtools
*/
class LazyTagStorage : TagStorage {
    this(ubyte[] chunk) {
        _chunk = chunk;
        if (std.system.endian != Endian.littleEndian) {
            fixByteOrder();
        }
    }

    final override Value opIndex(string key) {
        assert(key.length == 2);
        if (_chunk.length < 4)
            throw new RangeError();
        
       size_t offset = 0;
       while (offset < _chunk.length) {
           if (_chunk[offset .. offset + 2] == key) {
               offset += 2;
               return readValue(offset);
           } else {
               offset += 2;
               skipValue(offset);
           }
       }
       throw new RangeError();
    }

    final override int opApply(int delegate(ref string k, ref Value v) dg) {
        size_t offset = 0;
        while (offset < _chunk.length) {
            auto key = cast(string)_chunk[offset .. offset + 2];
            offset += 2;
            auto val = readValue(offset);
            auto res = dg(key, val);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

private:
    ubyte[] _chunk;

    Value readValue(ref size_t offset) {

        string readValueArrayTypeHelper() {
            char[] cases;
            foreach (c2t; ArrayElementTagValueTypes) {
                cases ~= 
                "case '"~c2t.ch~"':".dup~
                "  auto begin = offset;"~
                "  auto end = offset + length * "~c2t.ValueType.stringof~".sizeof;"~
                "  offset = end;"~
                "  return Value(cast("~c2t.ValueType.stringof~"[])(_chunk[begin .. end]));";
            }
            return to!string("switch (elem_type) {" ~ cases ~
                   "  default: throw new UnknownTagTypeException(to!string(elem_type));"~
                   "}");
        }

        string readValuePrimitiveTypeHelper() {
            char[] cases;
            foreach (c2t; PrimitiveTagValueTypes) {
                cases ~= "case '"~c2t.ch~"':"~
                         "  ubyte* p = _chunk.ptr + offset;"~ 
                         "  auto value = *(cast("~c2t.ValueType.stringof~"*)p);"~
                         "  offset += value.sizeof;"~
                         "  return Value(value);".dup;
            }
            return to!string("switch (type) {" ~ cases ~
                   "  default: throw new UnknownTagTypeException(to!string(type));"~
                   "}");
        }

        char type = cast(char)_chunk[offset++];
        if (type == 'Z' || type == 'H') {
            auto begin = offset;
            while (_chunk[offset++] != 0) {}
            // return string with stripped '\0'
            auto v = Value(cast(string)_chunk[begin .. offset - 1]);
            if (type == 'H') {
                v.setHexadecimalFlag();
            }
            return v;
        } else if (type == 'B') {
            char elem_type = cast(char)_chunk[offset++];
            uint length = *(cast(uint*)(_chunk.ptr + offset));
            offset += uint.sizeof;
            mixin(readValueArrayTypeHelper());
        } else {
            mixin(readValuePrimitiveTypeHelper());
        }
    }

    void skipValue(ref size_t offset) {
        char type = cast(char)_chunk[offset++];
        if (type == 'Z' || type == 'H') {
            while (_chunk[offset++] != 0) {}
        } else if (type == 'B') {
            char elem_type = cast(char)_chunk[offset++];
            auto length = *(cast(uint*)(_chunk.ptr + offset));
            offset += uint.sizeof + charToSizeof(elem_type) * length;
        } else {
            offset += charToSizeof(type);
        }
    }

    void fixByteOrder() {
        /* TODO: TEST ON BIG-ENDIAN SYSTEM!!! */
        ubyte* p = _chunk.ptr;
        ubyte* end = _chunk.ptr + _chunk.length;
        while (p < end) {
            p += 2; // skip tag name
            char type = *(cast(char*)p);
            ++p; // skip type
            if (type == 'Z' || type == 'H') {
                while (*p != 0) { // zero-terminated
                    ++p;          // string
                }
                ++p; // skip '\0'
            } else if (type == 'B') { // array
                char elem_type = *(cast(char*)p);
                uint size = charToSizeof(elem_type);
                switchEndianness(p, uint.sizeof);
                uint length = *(cast(uint*)p);
                p += uint.sizeof; // skip length
                if (size != 1) {
                    for (auto j = 0; j < length; j++) {
                        switchEndianness(p, size);
                        p += size;
                    }
                } else {
                    // skip 
                    p += length;
                }
            } else {
                uint size = charToSizeof(type);
                if (size != 1) {
                    switchEndianness(p, size);
                    p += size;
                } else {
                    ++p;
                }
            }
        }
    }
}
