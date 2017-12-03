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

struct RefSequence {
  size_d length;
  string name;
}

struct Header {
  string id;
  string text;
  RefSequence[] refs;

  this(this) {
    throw new Exception("Header does not have copy semantics");
  }

}

enum Offset {
  bin_mq_nl, flag_nc, l_seq, next_refID, next_pos, tlen, read_name, cigar, seq, qual, tag
};

/**
   Raw Read buffer containing unparsed data. It is read-only. When
   using fields beyond refid,pos use ProcessRead2 instead because it
   caches results.
*/

struct Read2 {
  uint refid;
  size_d pos;
  private ubyte[] _data;
  private uint[11] offset_table; // pre-computed offsets

  this(uint __refid, size_d __pos, ubyte[] data) {
    refid = __refid; pos = __pos; _data = data;
    offset_table = [Offset.bin_mq_nl*int.sizeof,Offset.flag_nc*int.sizeof,Offset.l_seq*int.sizeof,
                    Offset.next_refID*int.sizeof,Offset.next_pos*int.sizeof,Offset.tlen*int.sizeof,
                    Offset.read_name*int.sizeof,uint.max,uint.max,uint.max,uint.max];
    offset_table[Offset.cigar] =
      offset_table[Offset.read_name] + cast(uint)(_l_read_name * char.sizeof);
    offset_table[Offset.seq] =
      offset_table[Offset.cigar] + cast(uint)(_n_cigar_op * uint.sizeof);
    offset_table[Offset.qual] =
      offset_table[Offset.seq] + (_l_seq + 1)/2;
    offset_table[Offset.tag] = _qual_offset + _l_seq;
  }

  this(this) {
    throw new Exception("Read2 does not have copy semantics");
  }

  // uses raw offset
  @property @trusted nothrow private const T fetch_raw(T)(uint raw_offset) {
    ubyte[] buf = cast(ubyte[])_data[raw_offset..raw_offset+T.sizeof];
    return cast(const(T))buf.read!(T,Endian.littleEndian)();
  }

  // uses table offset
  @property @trusted nothrow private const T fetch(T)(uint offset) {
    auto raw_offset = offset_table[offset];
    return fetch_raw!T(raw_offset);
  }

  @property @trusted nothrow private const
  uint _bin_mq_nl()        { return fetch!uint(Offset.bin_mq_nl); }
  @property @trusted nothrow private const
  uint _flag_nc()          { return fetch!uint(Offset.flag_nc); }
  @property @trusted nothrow private const
  int sequence_length()    { return fetch!int(Offset.l_seq); }
  @property @trusted nothrow private const
  int _next_refID()        { return fetch!int(Offset.next_refID); }
  @property @trusted nothrow private const
  int _next_pos()          { return fetch!int(Offset.next_pos); }
  @property @trusted nothrow private const
  int _tlen()              { return fetch!int(Offset.tlen); }
  @property @trusted nothrow private const
  ushort _bin()            { return _bin_mq_nl >> 16; }
  @property @trusted nothrow private const
  ubyte _mapq()            { return (_bin_mq_nl >> 8) & 0xFF; }
  @property @trusted nothrow private const
  ubyte _l_read_name()     { return _bin_mq_nl & 0xFF; }
  @property @trusted nothrow private const
  ushort _flag()           { return _flag_nc >> 16; }
  @property @trusted nothrow private const
  ushort _n_cigar_op()     { return _flag_nc & 0xFFFF; }
  @property @trusted nothrow private const
  uint _read_name_offset() { return offset_table[Offset.read_name]; }
  @property @trusted nothrow private const
  uint _seq_offset()       { return offset_table[Offset.seq] ; }
  @property @trusted nothrow private const
  uint _qual_offset()      { return offset_table[Offset.qual] ; }
  @property @trusted nothrow private const
  uint _tags_offset()      { return offset_table[Offset.tag]; }
  @property @trusted nothrow private
  ubyte[] raw_sequence()   { return _data[offset_table[Offset.seq]..offset_table[Offset.qual]]; }

  alias sequence_length _l_seq;

  string toString() {
    return "<** " ~ Read2.stringof ~ " (data size " ~ to!string(_data.length) ~ ") " ~ to!string(refid) ~ ":" ~ to!string(pos) ~ " length " ~ to!string(sequence_length) ~ ">";
  }

}

/**
   ProcessRead2 provides a caching mechanism for Read2 fields. Use
   this when you need to access field/elements multiple times. Note
   that ProcessRead2 becomes invalid when Read2 goes out of scope.
*/
struct ProcessRead2 {
  private Read2 *_read2;
  Nullable!int sequence_length2;

  this(ref Read2 _r) {
    _read2 = cast(Read2 *)&_r;
  }

  @property @trusted int sequence_length() {
    if (sequence_length2.isNull)
      sequence_length2 = _read2.sequence_length;
    return sequence_length2;
  }

  @property ubyte[] raw_sequence() {
    return _read2.raw_sequence();
  }

  string toString() {
    return "<** " ~ ProcessRead2.stringof ~ ") " ~ to!string(_read2.refid) ~ ":" ~ to!string(_read2.pos) ~ " length " ~ to!string(sequence_length) ~ ">";
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
