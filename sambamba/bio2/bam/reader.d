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

struct RefSequence {
  size_t length;
  string name;
}

struct Header {
  string id;
  string text;
  RefSequence[] refs;

  this(this) {
    throw new Exception("struct has copy semantics");
  }

}

struct Read2 {
  uint refid;
  size_t pos;
  ubyte[] data;

  this(this) {
    throw new Exception("struct has copy semantics");
  }

  /*
  size_t pos() {
    data.read!(T,Endian.littleEndian)()
  }
  */

  string toString() {
    return "<**" ~ to!string(refid) ~ ":" ~ to!string(pos) ~ ">";
  }

}

struct BamReader2 {
  BgzfStream stream;
  Header header;

  this(string fn) {
    stream = BgzfStream(fn);
  }

  int opApply(scope int delegate(ref Read2) dg) {
    // parse the header
    ubyte[4] ubyte4;
    stream.read(ubyte4);
    enforce(ubyte4 == BAM_MAGIC,"Invalid file format: expected BAM magic number");
    immutable text_size = stream.read!int();
    immutable text = stream.read!string(text_size);
    header = Header(BAM_MAGIC,text);
    immutable n_refs = stream.read!int();
    foreach(int n_ref; 0..n_refs) {
      immutable l_name = stream.read!int();
      auto ref_name = stream.read!string(l_name);
      immutable l_ref = stream.read!int();
      header.refs ~= RefSequence(l_ref,ref_name);
    }
    // parse the reads
    while (!stream.eof()) {
      immutable block_size = stream.read!int();
      immutable refid = stream.read!int();
      immutable pos = stream.read!int();

      ubyte[] data = new ubyte[block_size-2*int.sizeof]; // Heap alloc
      auto read = Read2(refid,pos,stream.read(data));
      dg(read);
    }
    return 0;
  }
}
