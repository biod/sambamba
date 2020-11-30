/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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
/// BAM records may carry arbitrary information in tags.
/// $(BR)
/// $(D Value) type provides convenient way to work with this information.
///
/// Example:
/// --------------------------------
/// import bio.std.hts.bam.reader, bio.std.hts.bam.tagvalue;
/// ...
/// auto bam = new BamReader("file.bam");
/// Value v = bam.reads.front["MD"];
/// assert(v.is_string);
/// v = 5;
/// assert(v.is_signed);     // because 5 is of type int which is signed
/// assert(v == "5");        // converted to string and then compared
/// v = "abc";
/// assert(v.is_string);
/// v = [1, 2, 3];           // integer and float arrays are supported
/// assert(v.is_numeric_array);
/// v = [1.5f, 2.3f, 17.0f]; // double[] arrays must be converted to float[]
/// assert(v.is_numeric_array);
/// v = 5.6;
/// assert(v.is_float);
/// v = -17;
/// assert(v.is_signed);
/// ----------------------------------
module bio.std.hts.bam.tagvalue;

public import std.conv;
import std.typetuple;
import std.exception;
import std.format;
import std.array;
import bio.core.utils.format;

import bio.std.hts.thirdparty.msgpack;

struct CharToType(char c, T) {
    /** symbol */
    enum ch = c;

    /** type which corresponds to the symbol
        according to SAM/BAM specification 
    */
    alias T ValueType;    
}

/**
  Thrown in case of unrecognized tag type
 */
class UnknownTagTypeException : Exception {
    this(string msg) { super(msg); }
}


alias TypeTuple!(CharToType!('A', char),
                 CharToType!('c', byte),  
                 CharToType!('C', ubyte),
                 CharToType!('s', short), 
                 CharToType!('S', ushort),
                 CharToType!('i', int),   
                 CharToType!('I', uint),
                 CharToType!('f', float))       PrimitiveTagValueTypes;

alias TypeTuple!(CharToType!('Z', string), 
                 CharToType!('H', string))      StringTagValueTypes;

alias TypeTuple!(CharToType!('c', byte),  
                 CharToType!('C', ubyte),
                 CharToType!('s', short), 
                 CharToType!('S', ushort),
                 CharToType!('i', int),   
                 CharToType!('I', uint),
                 CharToType!('f', float))       ArrayElementTagValueTypes;

/*
  Useful in TagStorage implementations, for skipping elements

  Params:
    c =         primitive type identifier

  Returns: size of corresponding type in bytes
*/
uint charToSizeof(char c) {
    string charToSizeofHelper() {
        char[] cases;
        foreach (c2t; PrimitiveTagValueTypes) {
            cases ~= "case '"~c2t.ch~"':"~
                     "  return "~to!string(c2t.ValueType.sizeof)~";".dup;
        }
        return "switch (c) { " ~ cases.idup ~
               "  default: " ~
               "    throw new UnknownTagTypeException(to!string(c));"~ 
               "}";
    }
    mixin(charToSizeofHelper());
}

/*
  Pair of type and its ubyte identifier. 

  (Currently, ubyte is enough, but that might change in the future.)
*/
struct TypeId(T, ubyte id) {
    enum Id = id;
    alias T Type;
}

/*
  Structure of type identifier:

                              0                                   1   

                             primitive                          array/string
                 something         null/nothing             numeric         string
            numeric      char           0                   0              Z       H   
    integer        float                0              [see left           0       0 
unsigned   signed       0               0               branch]            0       0
 [ size in bytes]  [size in bytes]      0            [element size]        1       1

     (TypeId >> 5) == elementType.sizeof

*/
alias TypeTuple!(TypeId!(char,     0b001_00_1_00),
        
                 TypeId!(ubyte,    0b001_0_0000), 
                 TypeId!(ushort,   0b010_0_0000), 
                 TypeId!(uint,     0b100_0__0__0__0__0), 
/* Let's take                         4  u  i  n  s  p                  
   uint as an                            n  n  u  o  r                  
   example                            b  s  t  m  m  i                  
                                      y  i  e  e  e  m
                                      t  g  g  r  t  i
                                      e  n  e  i  h  t
                                      s  e  r  c  i  i
                                         d        n  v
                                                  g  e
*/   
 

                 TypeId!(byte,     0b001_1_0000),
                 TypeId!(short,    0b010_1_0000), 
                 TypeId!(int,      0b100_1_0000), 

                 TypeId!(float,    0b100_01_000),

                 TypeId!(ubyte[],  0b001_000_01),
                 TypeId!(ushort[], 0b010_000_01),
                 TypeId!(uint[],   0b100_000_01),

                 TypeId!(byte[],   0b001_100_01),
                 TypeId!(short[],  0b010_100_01),
                 TypeId!(int[],    0b100_100_01),

                 TypeId!(float[],  0b100_01_001),

                 TypeId!(string,   0b001_00_101),
                 TypeId!(string,   0b001_01_101),
                 TypeId!(typeof(null), 0b0000_0010))
    TypeIdMap;

private immutable hexStringTag = 0b001_01_101;

private template GetType(U) {
    alias U.Type GetType;
}

/// Get tag for type T.
///
/// Useful for comparison with tag field of Value struct.
/// 
/// Example:
/// -----------------------------------
/// Value v = "zzz";
/// assert(v.tag == GetTypeId!string);
/// -----------------------------------
template GetTypeId(T) {
    ///
    enum GetTypeId = TypeIdMap[staticIndexOf!(T, staticMap!(GetType, TypeIdMap))].Id;
}

string generateUnion() {
    char[] u = "union U {".dup;
    foreach (t; PrimitiveTagValueTypes) {
        u ~= t.ValueType.stringof ~ " " ~ t.ch ~ ";".dup;
    }
    foreach (t; StringTagValueTypes) {
        u ~= t.ValueType.stringof ~ " " ~ t.ch ~ ";".dup;
    }
    foreach (t; ArrayElementTagValueTypes) {
        u ~= t.ValueType.stringof ~ "[] " ~ 'B' ~ t.ch ~ ";".dup;
    }
    u ~= "}; U u;".dup;
    return u.idup;
}

template ArrayOf(T) {
    alias T[] ArrayOf;
}

string injectOpAssign() {
    char[] cs;

    foreach (t; PrimitiveTagValueTypes) {
        cs ~= "final void opAssign(" ~ t.ValueType.stringof ~ " value) {" ~
              "  this.u." ~ t.ch ~ " = value;" ~
              "  this._tag = " ~ to!string(GetTypeId!(t.ValueType)) ~ ";" ~
              "  this.bam_typeid = '" ~ t.ch ~ "';" ~
              "}";
    }

    cs ~= "final void opAssign(string value) {" ~
          "  this.u.Z = value;" ~
          "  this._tag = " ~ to!string(GetTypeId!string) ~ ";" ~
          "  this.bam_typeid = 'Z';" ~
          "}";

    foreach (t; ArrayElementTagValueTypes) {
        cs ~= "final void opAssign(" ~ t.ValueType.stringof ~ "[] value) {" ~
              "  this.u.B" ~ t.ch ~ " = value;" ~
              "  this._tag = " ~ to!string(GetTypeId!(ArrayOf!(t.ValueType))) ~ ";" ~
              "  this.bam_typeid = '" ~ t.ch ~ "';" ~
              "}";
    }

    return cs.idup;
}

string injectOpCast() {
    char[] cs = "static if".dup;

    string injectSwitchPrimitive(string requested_type) 
    {
        char[] cs = `switch (_tag) {`.dup;
              
        foreach (t2; PrimitiveTagValueTypes) {
            cs ~= `case GetTypeId!`~t2.ValueType.stringof~`: `~
                  `    return to!T(u.`~t2.ch~`);`.dup;
        }

        cs ~= `    default: throw new ConvException("Cannot convert Value to `~
                                                     requested_type~`");`~
              `}`;
        return cs.idup;
    }

    string injectSwitchArrayElement(string requested_type) 
    {
        char[] cs = `switch (_tag) {`.dup;
              
        foreach (t2; ArrayElementTagValueTypes) {
            cs ~= `case GetTypeId!(`~t2.ValueType.stringof~`[]): `~
                  `    return to!T(u.B`~t2.ch~`);`.dup;
        }

        cs ~= `    default: throw new ConvException("Cannot convert Value to `~
                                                     requested_type~`");`~
              `}`;
        return cs.idup;
    }

    foreach (t; TypeTuple!(byte, ubyte, short, ushort, int, uint,
                           char, float, double, real, long, ulong))
    {
        cs ~= `(is(T == `~t.stringof~`)) {`~
              injectSwitchPrimitive(t.stringof)~
              `} else static if`.dup;
    }

    foreach (t; ArrayElementTagValueTypes) {
        cs ~= `(is(T == ` ~ t.ValueType.stringof ~ `[])) {` ~
              injectSwitchArrayElement(t.ValueType.stringof ~ "[]")~
              `} else static if `;
    }

    cs ~= `(is(T == string)) {` ~
          `  if (is_string) {` ~
          `    return bam_typeid == 'Z' ? u.Z : u.H;`~
          `  } else if (is_integer || is_float || is_character) {`~
          `    `~injectSwitchPrimitive("string")~
          `  } else {`~
                 injectSwitchArrayElement("string")~
          `  }`~
          `}`.dup;

    return "final T opCast(T)() const {" ~ cs.idup ~ "}";
}

/**
  Struct for representing tag values. 

  Tagged union, allows to store 
  8/16/32-bit integers, floats, chars, strings, 
  and arrays of integers/floats.
*/
struct Value {

    /*
      Notice that having union first allows to do simple casts,
      without using opCast(). That's a bit hackish but
      allows for better speed.
     */
    private mixin(generateUnion());

    /**
      If this is an array, one of [cCsSiIf].
      Otherwise, one of [AcCsSiIfZH]

      See SAM/BAM specification for details.
    */
    public char bam_typeid;

    /*
                                    WARNING:

    Currently, type identifier for (u)int requires 8 bits.
    Fortunately, SAM/BAM specification doesn't use bigger integer types.
    However, in case of need to extend the hierarchy, the type
    should be changed from ubyte to something bigger. 
    */
    ubyte _tag;

    /// Designates the type of currently stored value.
    ///
    /// Supposed to be used externally for checking type with GetTypeId.
    ubyte tag() @property const {
        return _tag;
    }

    mixin(injectOpAssign());
    mixin(injectOpCast());

    ///
    final void opAssign(Value v) {
        bam_typeid = v.bam_typeid;
        _tag = v._tag;
        u = v.u;
    }

    /// ditto
    final void opAssign(typeof(null) n) {
        _tag = GetTypeId!(typeof(null));
    }

    ///
    final bool opEquals(T)(const T val) {
        try {
            return to!T(this) == val;
        } catch (ConvException e) {
            return false;
        }
    }

    ///
    string toString() const {
        return opCast!string();
    }

    ///
    this(T)(T value) {
        opAssign(value);
    }
 
    /// sets 'H' tag instead of default 'Z'. Is not expected to be used much.
    void setHexadecimalFlag() {

        enforce(this.is_string);
      
        bam_typeid = 'H';
        _tag = hexStringTag;

        if (_tag != 0b111) { 
            u.H = u.Z;
        }
    }

    /// Holds $(D null). Represents non-existing tag. Such values are used to remove tags.
    bool is_nothing() @property const { return _tag == GetTypeId!(typeof(null)); }

    /// char
    bool is_character() @property const { return _tag == GetTypeId!char; }

    /// float
    bool is_float() @property const { return _tag == GetTypeId!float; }

    /// ubyte[]/byte[]/ushort[]/short[]/uint[]/int[]/float[]
    bool is_numeric_array() @property const { return (_tag & 0b111) == 0b001; }

    /// ubyte[]/byte[]/ushort[]/short[]/uint[]/int[]
    bool is_array_of_integers() @property const { return (_tag & 0b1111) == 0b0001; }

    /// float[]
    bool is_array_of_floats() @property const { return (_tag & 0b1111) == 0b1001; }

    /// ubyte/byte/ushort/short/uint/int
    bool is_integer() @property const { return (_tag & 0b1111) == 0; }

    /// ubyte/ushort/uint
    bool is_unsigned() @property const { return (_tag & 0b11111) == 0; }

    /// byte/short/int
    bool is_signed() @property const { return (_tag & 0b11111) == 0b10000; }

    /// 'Z' or 'H' tag
    bool is_string() @property const { return (_tag & 0b111) == 0b101; }

    /// 'H' tag
    bool is_hexadecimal_string() @property const { return (_tag & 0b1101) == 0b1101; }

    /// Serializes value in MessagePack format
    public void toMsgpack(Packer)(ref Packer packer) const {
        switch (_tag) {
            case GetTypeId!byte: packer.pack(*cast(byte*)(&u)); break;
            case GetTypeId!ubyte: packer.pack(*cast(ubyte*)(&u)); break;
            case GetTypeId!short: packer.pack(*cast(short*)(&u)); break;
            case GetTypeId!ushort: packer.pack(*cast(ushort*)(&u)); break;
            case GetTypeId!int: packer.pack(*cast(int*)(&u)); break;
            case GetTypeId!uint: packer.pack(*cast(uint*)(&u)); break;

            case GetTypeId!float: packer.pack(*cast(float*)(&u)); break;
            case GetTypeId!string: packer.pack(*cast(char[]*)(&u)); break;
            case hexStringTag: packer.pack(*cast(char[]*)(&u)); break;
            case GetTypeId!char: packer.pack(*cast(ubyte*)(&u)); break;

            case GetTypeId!(byte[]): packer.pack(*cast(byte[]*)(&u)); break;
            case GetTypeId!(ubyte[]): packer.pack(*cast(ubyte[]*)(&u)); break;
            case GetTypeId!(short[]): packer.pack(*cast(short[]*)(&u)); break;
            case GetTypeId!(ushort[]): packer.pack(*cast(ushort[]*)(&u)); break;
            case GetTypeId!(int[]): packer.pack(*cast(int[]*)(&u)); break;
            case GetTypeId!(uint[]): packer.pack(*cast(uint[]*)(&u)); break;
            case GetTypeId!(float[]): packer.pack(*cast(float[]*)(&u)); break;

            case GetTypeId!(typeof(null)): packer.pack(null); break;
            default: break;
        }
    }

    /// SAM representation
    string toSam()() const {
        auto w = appender!(char[])();
        toSam((const(char)[] s) { w.put(s); });
        return cast(string)w.data;
    }

    /// ditto
    void toSam(Sink)(auto ref Sink sink) const 
        if (isSomeSink!Sink)
    {
        if (is_integer) {
            sink.write("i:");
            switch (_tag) {
                case GetTypeId!byte: sink.write(*cast(byte*)(&u)); break;
                case GetTypeId!ubyte: sink.write(*cast(ubyte*)(&u)); break;
                case GetTypeId!short: sink.write(*cast(short*)(&u)); break;
                case GetTypeId!ushort: sink.write(*cast(ushort*)(&u)); break;
                case GetTypeId!int: sink.write(*cast(int*)(&u)); break;
                case GetTypeId!uint: sink.write(*cast(uint*)(&u)); break;
                default: break;
            }
        } else if (is_numeric_array) {
            sink.write("B:");
            sink.write(bam_typeid);
            sink.write(',');
            switch (_tag) {
                case GetTypeId!(byte[]): sink.writeArray(*cast(byte[]*)(&u), ','); break;
                case GetTypeId!(ubyte[]): sink.writeArray(*cast(ubyte[]*)(&u), ','); break;
                case GetTypeId!(short[]): sink.writeArray(*cast(short[]*)(&u), ','); break;
                case GetTypeId!(ushort[]): sink.writeArray(*cast(ushort[]*)(&u), ','); break;
                case GetTypeId!(int[]): sink.writeArray(*cast(int[]*)(&u), ','); break;
                case GetTypeId!(uint[]): sink.writeArray(*cast(uint[]*)(&u), ','); break;
                case GetTypeId!(float[]): sink.writeArray(*cast(float[]*)(&u), ','); break;
                default: break;
            }
        } else {
            switch (_tag) {
                case GetTypeId!float: sink.write("f:"); sink.write(*cast(float*)(&u)); break;
                case GetTypeId!string: sink.write("Z:"); sink.write(*cast(const(char)[]*)(&u)); break;
                case hexStringTag: sink.write("H:"); sink.write(*cast(const(char)[]*)(&u)); break;
                case GetTypeId!char: sink.write("A:"); sink.write(*cast(char*)(&u)); break;
                default: break;
            }
        }
    }

    /// JSON representation
    string toJson()() const {
        auto w = appender!(char[])();
        toJson((const(char)[] s) { w.put(s); });
        return cast(string)w.data;
    }

    /// ditto
    void toJson(Sink)(auto ref Sink sink) const 
        if (isSomeSink!Sink)
    {
        switch (_tag) {
            case GetTypeId!byte: sink.writeJson(*cast(byte*)(&u)); break;
            case GetTypeId!ubyte: sink.writeJson(*cast(ubyte*)(&u)); break;
            case GetTypeId!short: sink.writeJson(*cast(short*)(&u)); break;
            case GetTypeId!ushort: sink.writeJson(*cast(ushort*)(&u)); break;
            case GetTypeId!int: sink.writeJson(*cast(int*)(&u)); break;
            case GetTypeId!uint: sink.writeJson(*cast(uint*)(&u)); break;
            case GetTypeId!(byte[]): sink.writeJson(*cast(byte[]*)(&u)); break;
            case GetTypeId!(ubyte[]): sink.writeJson(*cast(ubyte[]*)(&u)); break;
            case GetTypeId!(short[]): sink.writeJson(*cast(short[]*)(&u)); break;
            case GetTypeId!(ushort[]): sink.writeJson(*cast(ushort[]*)(&u)); break;
            case GetTypeId!(int[]): sink.writeJson(*cast(int[]*)(&u)); break;
            case GetTypeId!(uint[]): sink.writeJson(*cast(uint[]*)(&u)); break;
            case GetTypeId!(float[]): sink.writeJson(*cast(float[]*)(&u)); break;
            case GetTypeId!float: sink.writeJson(*cast(float*)(&u)); break;
            case GetTypeId!string: sink.writeJson(*cast(string*)(&u)); break;
            case hexStringTag: sink.writeJson(*cast(string*)(&u)); break;
            case GetTypeId!char: sink.writeJson(*cast(char*)(&u)); break;
            default: break;
        }
    }
}

Value readValueFromArray(char type, const(ubyte)[] bytes, ref size_t offset) {
    string readValueArrayTypeHelper() {
        char[] cases;
        foreach (c2t; ArrayElementTagValueTypes) {
            cases ~=
            "case '"~c2t.ch~"':".dup~
            "  auto begin = offset;"~
            "  auto end = offset + length * "~c2t.ValueType.stringof~".sizeof;"~
            "  offset = end;"~
            "  return Value(cast("~c2t.ValueType.stringof~"[])(bytes[begin .. end]));";
        }
        return to!string("switch (elem_type) {" ~ cases ~
               "  default: throw new UnknownTagTypeException(to!string(elem_type));"~
               "}");
    }

    string readValuePrimitiveTypeHelper() {
        char[] cases;
        foreach (c2t; PrimitiveTagValueTypes) {
            cases ~= "case '"~c2t.ch~"':"~
                     "  auto p = bytes.ptr + offset;"~
                     "  auto value = *(cast("~c2t.ValueType.stringof~"*)p);"~
                     "  offset += value.sizeof;"~
                     "  return Value(value);".dup;
        }
        return to!string("switch (type) {" ~ cases ~
               "  default: throw new UnknownTagTypeException(to!string(type));"~
               "}");
    }

    if (type == 'Z' || type == 'H') {
        auto begin = offset;
        while (bytes[offset++] != 0) {}
        // return string with stripped '\0'
        auto v = Value(cast(string)bytes[begin .. offset - 1]);
        if (type == 'H') {
            v.setHexadecimalFlag();
        }
        return v;
    } else if (type == 'B') {
        char elem_type = cast(char)bytes[offset++];
        uint length = *(cast(uint*)(bytes.ptr + offset));
        offset += uint.sizeof;
        mixin(readValueArrayTypeHelper());
    } else {
        mixin(readValuePrimitiveTypeHelper());
    }
}
