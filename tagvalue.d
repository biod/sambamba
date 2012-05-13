module tagvalue;

import std.variant;

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

import std.typetuple;

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
