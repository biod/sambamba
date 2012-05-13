module tagvalue;

import std.variant;
import std.typetuple;
import std.conv : to;

/**
  Represents tag value.

  Algebraic type, includes all types mentioned in SAM/BAM specification.

  */
alias Algebraic!(char, byte, ubyte, short, ushort, int, uint, float,
                 string, byte[], ubyte[], short[], ushort[],
                 int[], uint[], float[]) Value;
/*    A -> char
      c -> byte
      C -> ubyte
      s -> short
      S -> ushort
      i -> int
      I -> uint
      f -> float
      H -> string
      Z -> string
      B -> array of [cCsSiIf] */

/**
  Struct to be used in classes which implement TagStorage interface.
*/
struct CharToType(char c, T) {
    /** symbol */
    enum ch = c;

    /** type which corresponds to the symbol
        according to SAM/BAM specification 
    */
    alias T ValueType;    
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
  Thrown in case of unrecognized tag type
 */
class UnknownTagTypeException : Exception {
    this(string msg) { super(msg); }
}

/**
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
        return to!string("switch (c) { " ~ cases ~
               "  default: " ~
               "    throw new UnknownTagTypeException(to!string(c));"~ 
               "}");
    }
    mixin(charToSizeofHelper());
}
