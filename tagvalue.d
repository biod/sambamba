module tagvalue;

import std.conv : to;
import std.typetuple;

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
    Currently, ubyte is more than enough 
    (there're only 17 tag types at the moment)
*/
struct TypeId(T, ubyte id) {
    enum Id = id;
    alias T Type;
}

/*
  Structure of type identifier:

                      0                                  1   

                  primitive                       array/string
            numeric      char                numeric         string
    integer        float              (see left branch)     Z       H   
unsigned   signed              
 [ size in bytes]               

*/
alias TypeTuple!(TypeId!(char,     0b00000_1_0),
        
                 TypeId!(ubyte,    0b001_0_000), 
                 TypeId!(ushort,   0b010_0_000), 
                 TypeId!(uint,     0b100_0_000), 

                 TypeId!(byte,     0b001_1_000),
                 TypeId!(short,    0b010_1_000), 
                 TypeId!(int,      0b100_1_000), 

                 TypeId!(float,    0b0000_1_00),

                 TypeId!(ubyte[],  0b001_00_01),
                 TypeId!(ushort[], 0b010_00_01),
                 TypeId!(uint[],   0b100_00_01),

                 TypeId!(byte[],   0b001_10_01),
                 TypeId!(short[],  0b010_10_01),
                 TypeId!(int[],    0b100_10_01),

                 TypeId!(float[],  0b0000_1_01),

                 TypeId!(string,   0b0000_0_11),
                 TypeId!(string,   0b0000_1_11)) 
    TypeIdMap;


template GetTypeId(T) {
    template GetType(U) {
        alias U.Type GetType;
    }
    private enum index = staticIndexOf!(T, staticMap!(GetType, TypeIdMap));
    enum GetTypeId = TypeIdMap[index].Id;
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
    return to!string(u);
}

template ArrayOf(T) {
    alias T[] ArrayOf;
}

string defineOpAssign() {
    char[] cs;

    foreach (t; PrimitiveTagValueTypes) {
        cs ~= "final void opAssign(" ~ t.ValueType.stringof ~ " value) {" ~
              "  this.u." ~ t.ch ~ " = value;" ~
              "  this.tag = " ~ to!string(GetTypeId!(t.ValueType)) ~ ";" ~
              "  this.bam_typeid = '" ~ t.ch ~ "';" ~
              "}";
    }

    cs ~= "final void opAssign(string value) {" ~
          "  this.u.Z = value;" ~
          "  this.tag = " ~ to!string(GetTypeId!string) ~ ";" ~
          "  this.bam_typeid = 'Z';" ~
          "}";

    foreach (t; ArrayElementTagValueTypes) {
        cs ~= "final void opAssign(" ~ t.ValueType.stringof ~ "[] value) {" ~
              "  this.u.B" ~ t.ch ~ " = value;" ~
              "  this.tag = " ~ to!string(GetTypeId!(ArrayOf!(t.ValueType))) ~ ";" ~
              "  this.bam_typeid = '" ~ t.ch ~ "';" ~
              "}";
    }

    return to!string(cs);
}

struct Value {
    /**
      If this is an array, one of [cCsSiIf].
      Otherwise, one of [AcCsSiIfZH]

      See SAM/BAM specification for details.
    */
    public char bam_typeid;

    private ubyte tag; /// for internal use

    private mixin(generateUnion());

    mixin(defineOpAssign());

    final void opAssign(Value v) {
        bam_typeid = v.bam_typeid;
        tag = v.tag;
        u = v.u;
    }

    this(T)(T value) {
        opAssign(value);
    }
 
    /// sets 'H' tag instead of default 'Z'. Is not expected to be used much.
    void setHexadecimalFlag() {
        bam_typeid = 'H';
        tag = 0b111;
        u.H = u.Z;
    }

    bool is_char() @property { return tag == GetTypeId!char; }
    bool is_float() @property { return tag == GetTypeId!float; }
    bool is_numeric_array() @property { return (tag & 0b11) == 0b01; }
    bool is_array_of_integers() @property { return (tag & 0b111) == 0b001; }
    bool is_array_of_floats() @property { return (tag & 0b111) == 0b101; }
    bool is_integer() @property { return (tag & 0b111) == 0; }

    /// true if the value is unsigned integer
    bool is_unsigned() @property { return (tag & 0b1111) == 0; }

    /// true if the value is signed integer
    bool is_signed() @property { return (tag & 0b1111) == 0b1000; }

    /// true if the value represents 'Z' or 'H' tag
    bool is_string() @property { return (tag & 0b11) == 0b11; }

    /** representation in SAM format
     
        Example:
        ----------
        Value v = 2.7;
        assert(v.to_sam() == "f:2.7");

        v = [1, 2, 3];
        assert(v.to_sam() == "B:i,1,2,3");
    */
    string to_sam() @property {
        if (this.is_numeric_array) {
            string toSamNumericArrayHelper() {
                char[] cases;
                foreach (t; ArrayElementTagValueTypes) 
                    cases ~= `case '`~t.ch~`':` ~
                             `  char[] str = "B:`~t.ch~`".dup;`~
                             `  foreach (elem; this.u.B`~t.ch~`) {`~
                             `    str ~= ',' ~ to!string(elem);`~
                             `  }`~
                             `  return cast(string)str;`.dup;
                return "switch (bam_typeid) { " ~ cases.idup ~ "default: assert(0); }";
            }
            mixin(toSamNumericArrayHelper());
        }
        if (this.is_integer) {
            switch (bam_typeid) {
                /* TODO: this probably can be optimized 
                         using some bit hacks
                         since we can compute sizeof in O(1)
                */
                case 'c': return "i:" ~ to!string(u.c);
                case 'C': return "i:" ~ to!string(u.C);
                case 's': return "i:" ~ to!string(u.s);
                case 'S': return "i:" ~ to!string(u.S);
                case 'i': return "i:" ~ to!string(u.i);
                case 'I': return "i:" ~ to!string(u.I);
                default: assert(0);
            }
        }
        if (this.is_float) {
            return "f:" ~ to!string(u.f);
        }
        switch (bam_typeid) {
            case 'Z': return "Z:" ~ u.Z;
            case 'H': return "H:" ~ u.H;
            case 'A': return "A:" ~ u.A;
            default: assert(0);
        }
    }
}

unittest {
    import std.stdio;

    writeln("Testing converting tags to SAM representation...");
    Value v = 5;
    assert(v.is_integer);
    assert(v.to_sam == "i:5");
    v = "abc";
    assert(v.is_string);
    assert(v.to_sam == "Z:abc");
    v = [1, 2, 3];
    assert(v.is_numeric_array);
    assert(v.to_sam == "B:i,1,2,3");
    v = [1.5, 2.3, 17.0];
    assert(v.is_numeric_array);
    assert(v.to_sam == "B:f,1.5,2.3,17");
    v = 5.6;
    assert(v.is_float);
    assert(v.to_sam == "f:5.6");
    v = -17;
    assert(v.is_signed);
    assert(v.to_sam == "i:-17");
    v = 297u;
    assert(v.is_unsigned);
    assert(v.to_sam == "i:297");

    short[] array_of_shorts = [4, 5, 6];
    v = array_of_shorts;
    assert(v.is_numeric_array);
    assert(v.to_sam == "B:s,4,5,6");
}
