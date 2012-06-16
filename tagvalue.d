module tagvalue;

public import std.conv;
import std.typetuple;
import std.exception;

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

/**
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

/**
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

                 TypeId!(float,    0b100_0_1_000),

                 TypeId!(ubyte[],  0b001_000_01),
                 TypeId!(ushort[], 0b010_000_01),
                 TypeId!(uint[],   0b100_000_01),

                 TypeId!(byte[],   0b001_010_01),
                 TypeId!(short[],  0b010_010_01),
                 TypeId!(int[],    0b100_010_01),

                 TypeId!(float[],  0b100_00_1_01),

                 TypeId!(string,   0b001_00_0_11),
                 TypeId!(string,   0b001_00_1_11),
                 TypeId!(typeof(null), 0b0000_0010))
    TypeIdMap;


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
    char[] cs = `static if `.dup;
    
    foreach (t; PrimitiveTagValueTypes) {
        cs ~= `(is(T == `~t.ValueType.stringof~`)) {`~
              `  if (this._tag != `~to!string(GetTypeId!(t.ValueType))~`) {`~
              `    throw new ConvException("Cannot convert Value to `~
                                           t.ValueType.stringof~`");`~
              `  }`~
              `  return this.u.`~t.ch~`;`~
              `} else static if `;
    }

    foreach (t; ArrayElementTagValueTypes) {
        cs ~= `(is(T == ` ~ t.ValueType.stringof ~ `[])) {` ~
              `  if (this._tag != `~to!string(GetTypeId!(ArrayOf!(t.ValueType)))~`) {`~
              `    throw new ConvException("Cannot convert Value to `~
                                           t.ValueType.stringof~`[]");`~
              `  }`~
              `  return this.u.B`~t.ch~`;`~
              `} else static if `;
    }

    cs ~= `(is(T == string)) {` ~
          `  if (is_string) {`
          `    return bam_typeid == 'Z' ? u.Z : u.H;`~
          `  } else {`~
          `    throw new ConvException("Cannot convert Value to string");`~
          `  }`~
          `}`.dup;

    return "final T opCast(T)() {" ~ cs.idup ~ "}";
}

/**
  Struct for representing tag values. 

  Tagged union, allows to store 
  8/16/32-bit integers, floats, chars, strings, 
  and arrays of integers/floats.

  Currently, opCast is very restrictive and requires that 
  the requested type is exactly the same as stored in Value
  (otherwise, ConvException is thrown). That means that
  you can't cast Value to string when it contains integer,
  although it's possible to convert integer to string.
*/
struct Value {

    /**
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
    ubyte tag() @property {
        return _tag;
    }

    mixin(injectOpAssign());
    mixin(injectOpCast());

    final void opAssign(Value v) {
        bam_typeid = v.bam_typeid;
        _tag = v._tag;
        u = v.u;
    }

    final void opAssign(typeof(null) n) {
        _tag = GetTypeId!(typeof(null));
    }

    final bool opEqual(T)(T val) {
        return to!T(this) == val;
    }

    /// Conversion to string occurs only when Value stores 
    /// 'Z' or 'H' tag. Otherwise ConvException is thrown.
    string toString() {
        return opCast!string();
    }

    this(T)(T value) {
        opAssign(value);
    }
 
    /// sets 'H' tag instead of default 'Z'. Is not expected to be used much.
    void setHexadecimalFlag() {

        enforce(this.is_string);

        bam_typeid = 'H';
        _tag = 0b111;
        u.H = u.Z;
    }

    bool is_nothing() @property { return _tag == GetTypeId!(typeof(null)); }

    bool is_character() @property { return _tag == GetTypeId!char; }
    bool is_float() @property { return _tag == GetTypeId!float; }
    bool is_numeric_array() @property { return (_tag & 0b11) == 0b01; }
    bool is_array_of_integers() @property { return (_tag & 0b111) == 0b001; }
    bool is_array_of_floats() @property { return (_tag & 0b111) == 0b101; }
    bool is_integer() @property { return (_tag & 0b1111) == 0; }

    /// true if the value is unsigned integer
    bool is_unsigned() @property { return (_tag & 0b11111) == 0; }

    /// true if the value is signed integer
    bool is_signed() @property { return (_tag & 0b11111) == 0b10000; }

    /// true if the value represents 'Z' or 'H' tag
    bool is_string() @property { return (_tag & 0b11) == 0b11; }

    /// true if the value represents 'H' tag
    bool is_hexadecimal_string() @property { return (_tag & 0b111) == 0b111; }
}
