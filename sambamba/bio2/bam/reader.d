/**
   New style BAM reader
*/

module sambamba.bio2.bam.reader;

import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.exception;
import std.file;
import std.stdio;
import std.string;
import std.typecons;
import std.bitmanip;

import bio.bam.constants;

import sambamba.bio2.bgzf;

struct Header2 {
  this(immutable(ubyte[]) buf) {
  }
}

alias Nullable!Header2 Header;

struct Read2 {
  immutable(ubyte[]) buffer;

  this(immutable(ubyte[]) buf) {
    buffer = buf;
  }

  string toString() {
    return to!string(buffer[32..34]);
  }
}

struct BamReader2 {

  BgzfBlocks bgzfblocks;
  Header header;

  this(string fn) {
    bgzfblocks = BgzfBlocks(fn);
  }

  int opApply(scope int delegate(Read2) dg) {
    foreach(immutable(ubyte[]) block; bgzfblocks) {
      stderr.writeln("Reading new block with length ",block.length);
      auto nblock = block[0..$];
      if (header.isNull) {
        // Start parsing the header
        ubyte[4] ubyte4;
        enforce(block[0..4] == BAM_MAGIC,"Invalid file format: expected BAM magic number");
        nblock = block[4..$];
        immutable text_size = nblock.read!(int,Endian.littleEndian)();
        stderr.writeln(text_size);
        nblock = nblock[text_size..$];
        immutable n_refs = nblock.read!(int,Endian.littleEndian)();
        foreach(int n_ref; 0..n_refs) {
          immutable l_name =  nblock.read!(int,Endian.littleEndian)();
          nblock = nblock[l_name..$];
          immutable l_ref =  nblock.read!(int,Endian.littleEndian)();
        }

        auto header_size = nblock.ptr - block.ptr;
        stderr.writeln("Header size ",header_size," nblock ",nblock.length);
        header = Header2(block[0..header_size]);
      }
      else {
        immutable block_size = nblock.read!(int,Endian.littleEndian)();
        stderr.writeln("Read block size ",block_size);
        auto read = Read2(nblock[0..block_size]);
        dg(read);
      }
    }
    return 0;
  }

}
