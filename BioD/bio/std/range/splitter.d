/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.range.splitter;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.stdio;

import std.range.primitives;

immutable ubyte[] SPLIT_ON = [ 0x20, 0x09, 0x0A, ';', ',' ];

/**
   SimpleSplitConv takes a range R (typically a text line) and splits
   it/tokenizes it on a list of characters. Essentially fields/tokens
   are split by tabs, semi-colons or comma's and spaces. This compares
   to C's strtok(str, ", \t;").

   This routine happens often in bioinformatics and is a replacement
   for the much unsafer C strtok.  This edition should also handle
   UTF.

   The default is to split on space, newline, tab, semi-colon and
   comma.
*/

struct SimpleSplitConv(R)
  if (isInputRange!R)
{
  R list, split_on;

  this(R range, R splits_on = cast(R)SPLIT_ON) {
    list = range;
    split_on = splits_on;
  }

  int opApply(scope int delegate(R) dg) {
    size_t start = 0;
    bool in_whitespace = false;
    foreach(size_t pos, c; list) {
      if (canFind(split_on,c)) { // hit split char
        if (!in_whitespace) { // emit
          auto token = list[start..pos];
          dg(token);
        }
        start = pos+1;
        in_whitespace = true;
      } else {
        in_whitespace = false;
      }
    }
    if (!in_whitespace) { // emit final
      auto token = list[start..$];
      dg(token);
    }
    return 0;
  }
}

unittest {
  auto s = cast(ubyte[])"hello 1 2 \t3  4 \n";
  assert(array(SimpleSplitConv!(ubyte[])(s)) == ["hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"  hello, 1 2 \t3  4 \n")) == ["","hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"hello, 1 2 \n\t3  4 \n")) == ["hello","1","2","3","4"]);
  assert(array(SimpleSplitConv!(ubyte[])(cast(ubyte[])"chr1:55365,55365,1")) == ["chr1:55365","55365","1"]);
}

/*
   Dirty fast_splitter is 3x faster than above elegant version. It does no heap
   allocations.
*/
R[] fast_splitter(R)(R[] tokens, R range, R splits_on = cast(R)SPLIT_ON) @nogc {
  // R[] tokens = new R[range.length]; // pre-allocate optimistially
  auto j = 0, prev_j = 0;
  bool in_whitespace = false;
  auto token_num = 0;
  for (; j<range.length ;) {
    bool found = false;
    auto check = range[j];
    foreach (c ; splits_on) {
      if (c==check) {
        found = true;
        break;
      }
    }
    if (found) {
      if (!in_whitespace) {
        tokens[token_num] = range[prev_j..j];
        token_num++;
      }
      prev_j = j+1;
      in_whitespace = true;
    }
    else {
      in_whitespace = false;
    }
    j++;
  }
  if (!in_whitespace) { // emit final
    tokens[token_num] = range[prev_j..$];
    token_num++;
  }
  // tokens.length = token_num;
  return tokens[0..token_num];
}

/*
   Same as above, but with one single heap allocation - it may be slightly
   slower.
*/
R[] fast_splitter(R)(R range, R splits_on = cast(R)SPLIT_ON) {
  R[] tokens = new R[range.length];
  return fast_splitter(tokens,range,splits_on);
}

unittest {
  auto s = "hello 1 2 \t3  4 \n";
  string[16] tokens; // preset buffer
  assert(fast_splitter(tokens,s) == ["hello", "1", "2", "3", "4"]);
  assert(fast_splitter("  hello, 1 2 \t3  4 \n") == ["","hello","1","2","3","4"]);
  assert(fast_splitter("hello, 1 2 \n\t3  4 \n") == ["hello","1","2","3","4"]);
  assert(fast_splitter(tokens,"chr1:55365,55365,1") == ["chr1:55365","55365","1"]);
}
