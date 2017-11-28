/**
   New style BAM reader
*/

import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.exception;
import std.file;
import std.stdio;
import std.string;
import std.bitmanip;

import sambamba.bio2.bgzf;

struct Read2 {
}

struct BamReader2 {

  BgzfBlocks bgzfblocks;

  this(string fn) {
    bgzfblocks = BgzfBlocks(fn);
  }

  int opApply(int delegate(ref int) operations) const  {
    foreach(ubyte[] block; bgzfblocks) {
      dg(block);
    }
    return 0;
  }

}
