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
import sambamba.bio2.constants;

// alias ulong size_d;

struct RefSequence {
  size_d length;
  string name;
}

struct Header {
  string id;
  string text;
  RefSequence[] refs;

  this(this) {
    throw new Exception("Header has copy semantics");
  }

}

enum Offset { l_seq = 2*int.sizeof,
              next_refID = 3*int.sizeof
};

/**
   Raw Read buffer containing unparsed data. It is read-only. When
   using fields beyond refid,pos use ProcessRead2 instead because it
   caches results.
*/

struct Read2 {
  uint refid;
  size_d pos;
  ubyte[] data;
  int refcount = 0;

  this(this) {
    throw new Exception("Read2 has copy semantics");
  }

  nothrow @property @trusted const T fetch(T)(size_t offset) {
    // this may be a bit slower than the original, but we'll cache anyway
    ubyte[] buf = cast(ubyte[])data[offset..offset+T.sizeof];
    return cast(const(T))buf.read!(T,Endian.littleEndian)();
  }

  nothrow @property @trusted const int sequence_length() {
    return fetch!int(Offset.l_seq);
  }

  string toString() {
    return "<** " ~ Read2.stringof ~ " (data size " ~ to!string(data.length) ~ ") " ~ to!string(refid) ~ ":" ~ to!string(pos) ~ " length " ~ to!string(sequence_length) ~ ">";
  }

}

/**
   ProcessRead2 provides a caching mechanism for Read2 fields. Use
   this when you need to access field/elements multiple times. Note
   that ProcessRead2 becomes invalid when Read2 goes out of scope.
*/
struct ProcessRead2 {
  Read2 *read2;
  Nullable!int sequence_length2;

  this(ref Read2 _r) {
    read2 = cast(Read2 *)&_r;
  }

  nothrow @property @trusted int sequence_length() {
    if (sequence_length2.isNull)
      sequence_length2 = read2.sequence_length;
    return sequence_length2;
  }

  string toString() {
    return "<** " ~ ProcessRead2.stringof ~ ") " ~ to!string(read2.refid) ~ ":" ~ to!string(read2.pos) ~ " length " ~ to!string(sequence_length) ~ ">";
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
