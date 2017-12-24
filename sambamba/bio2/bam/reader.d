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

import bio.bam.cigar;
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

  @disable this(this); // disable copy semantics;

}

void fetch_bam_header(ref Header header, ref BgzfStream stream) {
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
}

template ReadFlags(alias flag) {
  @property bool is_paired()                nothrow { return cast(bool)(flag & 0x1); }
  /// Each segment properly aligned according to the aligner
  @property bool is_proper_pair()           nothrow { return cast(bool)(flag & 0x2); }
  @property bool is_unmapped()              nothrow { return cast(bool)(flag & 0x4); }
  @property bool is_mapped()                nothrow { return cast(bool)(!(flag & 0x4)); }
  @property bool mate_is_unmapped()         nothrow { return cast(bool)(flag & 0x8); }
  @property bool is_reverse_strand()        nothrow { return cast(bool)(flag & 0x10); }
  @property bool mate_is_reverse_strand()   nothrow { return cast(bool)(flag & 0x20); }
  @property bool is_first_of_pair()         nothrow { return cast(bool)(flag & 0x40); }
  @property bool is_second_of_pair()        nothrow { return cast(bool)(flag & 0x80); }
  @property bool is_secondary_alignment()   nothrow { return cast(bool)(flag & 0x100); }
  @property bool is_qc_fail() nothrow {
    assert(is_mapped());
    return cast(bool)(flag & 0x200); }
  alias is_qc_fail failed_quality_control;
  /// PCR or optical duplicate
  @property bool is_duplicate()             nothrow { return cast(bool)(flag & 0x400); }
  /// Supplementary alignment
  @property bool is_supplementary()         nothrow { return cast(bool)(flag & 0x800); }
}

template CheckMapped(alias refid) {
  @property nothrow bool is_unmapped2() {
    return is_unmapped;
  }
  @property nothrow bool is_mapped2() {
    debug {
      if (is_mapped) {
        assert(refid != -1, "ref_id can not be -1 for mapped read");  // BAM spec
      }
    }
    return !is_unmapped;
  }
}

enum Offset {
  bin_mq_nl=0, flag_nc=4, l_seq=8, next_refID=12, next_pos=16, tlen=20, read_name=24
};

/**
   Raw Read buffer containing unparsed data. It should be considered
   read-only.

   Current implementation is a cluct (class-struct hybrid). The _data
   pointer is shared when ReadBlob is assigned to another variable
   (i.e., there is a remote dependency). The advantage is that for
   each Read data gets allocated on the heap only once.

   All offsets are indexed on init (except for tags).  When using
   fields beyond refid,pos use ProcessReadBlob instead because it
   caches values.
*/

struct ReadBlob {
  RefId refid;   // -1 is invalid (BAM Spec)
  GenomePos pos; // 0 coordinate based (BAM spec)
  private ubyte[] _data;
  uint offset_cigar=int.max, offset_seq=int.max, offset_qual=int.max;

  mixin ReadFlags!(_flag_nc);
  mixin CheckMapped!(refid);

  /*
  this(RefId ref_id, GenomePos read_pos, ubyte[] buf) {
    refid = ref_id;
    pos = read_pos;
    _data = buf;
  }
  */

  // Turn ReadBlob into class-struct hybrid or a cluct ;)
  // @disable this(this); // disable copy semantics;

  @property @trusted nothrow private const T fetch(T)(uint raw_offset) {
    ubyte[] buf = cast(ubyte[])_data[raw_offset..raw_offset+T.sizeof];
    return cast(const(T))buf.read!(T,Endian.littleEndian)();
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
  int _tlen()              { return fetch!int(Offset.tlen); } // avoid using TLEN
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
  uint _read_name_offset() { return Offset.read_name; }
  @property @trusted nothrow private
  uint _cigar_offset()       {
    if (offset_cigar == int.max)
      offset_cigar = Offset.read_name + cast(uint)(_l_read_name * char.sizeof);
    return offset_cigar;
  }
  @property @trusted nothrow private
  uint _seq_offset()       {
    if (offset_seq == int.max)
      offset_seq = _cigar_offset + cast(uint)(_n_cigar_op * uint.sizeof);
    return offset_seq;
  }
  @property @trusted nothrow private
  uint _qual_offset()      {
    if (offset_qual == int.max)
      offset_qual = _seq_offset + (sequence_length + 1)/2;
    return offset_qual;
  }
  @property @trusted nothrow private
  uint _tags_offset()      { return _qual_offset + sequence_length; }
  @property @trusted nothrow private
  ubyte[] read_name()      { return _data[_read_name_offset.._cigar_offset]; }
  @property @trusted nothrow private
  ubyte[] raw_cigar()      { return _data[_cigar_offset.._seq_offset]; }
  @property @trusted nothrow private
  ubyte[] raw_sequence()   { return _data[_seq_offset.._qual_offset]; }

  alias sequence_length _l_seq;
  alias _mapq mapping_quality; // MAPQ

  string toString() {
    return "<** " ~ ReadBlob.stringof ~ " (data size " ~ to!string(_data.length) ~ ") " ~ to!string(refid) ~ ":" ~ to!string(pos) ~ " length " ~ to!string(sequence_length) ~ ">";
  }

}

/**
   ProcessReadBlob provides a caching mechanism for ReadBlob fields. Use
   this when you need to access field/elements multiple times. Note
   that ProcessReadBlob becomes invalid when ReadBlob goes out of scope.
*/
struct ProcessReadBlob {
  private Nullable!ReadBlob _read2;
  Nullable!int sequence_length2;
  private Nullable!string sequence2, read_name2;
  private Nullable!CigarOperations cigar2;
  private Nullable!GenomePos consumed_reference_bases2;

  mixin ReadFlags!(_flag);
  mixin CheckMapped!(ref_id);

  this(Nullable!ReadBlob _r) {
    _read2 = _r;
  }

  @property nothrow bool isNull() {
    return _read2.isNull;
  }

  @property nothrow RefId ref_id() {
    assert(_read2.is_mapped2,"Trying to get ref_id an unmapped read");
    return _read2.refid;
  }

  private @property nothrow RefId _flag() {
    return _read2._flag_nc;
  }

  alias ref_id refid;

  /// Get the start position on the reference sequence (better use start_loc)
  @property GenomePos start_pos() {
    assert(_read2.is_mapped2,"Trying to get pos on an unmapped read"); // BAM spec
    enforce(_read2.pos < GenomePos.max);
    return cast(GenomePos)_read2.pos;
  }

  /// Get the end position on the reference sequence (better use end_loc)
  @property GenomePos end_pos() {
    assert(sequence_length > 0, "Trying to get end_pos on an empty read sequence");
    assert(!consumed_reference_bases.isNull);
    return start_pos + consumed_reference_bases;
  }

  @property GenomeLocation start_loc() {
    return GenomeLocation(ref_id,start_pos);
  }

  @property GenomeLocation end_loc() {
    return GenomeLocation(ref_id,end_pos);
  }

  @property @trusted MappingQuality mapping_quality() { // MAPQ
    assert(_read2.is_mapped2,"Trying to get MAPQ on an unmapped read"); // BAM spec
    return MappingQuality(_read2.mapping_quality);
  }

  @property @trusted int tlen() { // do not use
    return _read2._tlen;
  }

  @property @trusted GenomePos sequence_length() {
    if (sequence_length2.isNull)
      sequence_length2 = _read2.sequence_length;
    return sequence_length2;
  }

  @property @trusted Nullable!GenomePos consumed_reference_bases() {
    if (consumed_reference_bases2.isNull) {
      assert(_read2.is_mapped2,"Trying to get consumed bases on an unmapped read"); // BAM spec
      assert(!read_name.isNull,"Trying to get CIGAR on RNAME is '*'"); // BAM spec
      auto raw = cast(uint[]) _read2.raw_cigar();
      if (raw.length==1 && raw[0] == '*')
        return consumed_reference_bases2; // null
      else {
        GenomePos bases = 0;
        for (size_t i = 0; i < raw.length; i++) {
          auto cigarop = CigarOperation(raw[i]);
          if (cigarop.is_query_consuming)
            bases += cigarop.length;
        }
        consumed_reference_bases2 = bases;
      }
    }
    return consumed_reference_bases2;
  }

  /// Return read name as a string. If unavailable returns
  /// null. Caches name.
  @property Nullable!string read_name() {
    if (read_name2.isNull) {
      assert(_read2.is_mapped2,"Trying to get RNAME on an unmapped read"); // BAM spec
      auto raw = _read2.read_name;
      if (raw.length == 0 || (raw.length ==1 && raw[0] == '*'))
        return read_name2; // null
      assert(raw.length < 255); // BAM spec
      if (raw[raw.length-1] == 0) // strip trailing C zero
        raw.length -= 1;
      read_name2 = Nullable!string(cast(string)raw);
    }
    return read_name2;
  }

  /// Returns Cigar as an array of operations. Returns null if no
  /// operations. Caches Cigar when there are operations.
  @property Nullable!CigarOperations cigar() {
    if (cigar2.isNull) {
      assert(_read2.is_mapped2,"Trying to get CIGAR on an unmapped read"); // BAM spec
      assert(!read_name.isNull,"Trying to get CIGAR on RNAME is '*'"); // BAM spec
      auto raw = cast(uint[]) _read2.raw_cigar();
      if (raw.length==0 || (raw.length==1 && raw[0] == '*'))
        return cigar2; // null
      else {
        auto s = new CigarOperation[raw.length]; // Heap alloc
        s.length = 0;
        for (size_t i = 0; i < raw.length; i++) {
          s ~= CigarOperation(raw[i]);
        }
        cigar2 = s;
      }
    }
    return cigar2;
  }

  /// Return human readable sequence fragment - null if
  /// undefined. Caches sequence.
  @property Nullable!string sequence() {
    if (sequence2.isNull) { // is it cached in sequence2?
      auto raw = _read2.raw_sequence();
      if (raw[0] == '*') {
        assert(raw.length == 1);
        return sequence2; // null
      }
      auto raw_length = (sequence_length + 1) / 2;
      char[16] convert = "=ACMGRSVTWYHKDBN";
      string s;
      s.reserve(sequence_length); // Heap alloc
      for (size_t i = 0; i < sequence_length; i++) {
        auto is_odd = i % 2;
        auto nuc = (is_odd ? raw[i/2] & 0b00001111 : (raw[i/2] & 0b11110000) >> 4);
        s ~= convert[nuc];
      }
      sequence2 = s;
    }
    return sequence2;
  }

  string toString() {
    return "<** " ~ ProcessReadBlob.stringof ~ ") " ~ to!string(_read2.refid) ~ ":" ~ to!string(_read2.pos) ~ " length " ~ to!string(sequence_length) ~ ">";
  }

}

/**
   BamReader2 is used for foreach loops
*/

struct BamReadBlobs {
  BgzfStream stream;
  Header header;

  this(string fn) {
    stream = BgzfStream(fn);
  }

  int opApply(scope int delegate(ref ReadBlob) dg) {
    fetch_bam_header(header, stream);
    // parse the reads
    while (!stream.eof()) {
      immutable block_size = stream.read!int();
      immutable refid = stream.read!int();
      immutable pos = stream.read!int();

      ubyte[] data = new ubyte[block_size-2*int.sizeof]; // Heap alloc FIXME
      auto read = ReadBlob(refid,pos,stream.read(data));
      dg(read);
    }
    return 0;
  }
}

/**
   Read streamer - use on single thread only
*/

// import core.memory : pureMalloc;

struct BamReadBlobStream {
  BgzfStream stream;
  Header header;
  Nullable!ReadBlob current;
  ubyte[] data; // in sync with current

  this(string fn) {
    stream = BgzfStream(fn);
    fetch_bam_header(header, stream);
    popFront();
  }

  bool empty() @property {
    return stream.eof();
  }

  ReadBlob front() {
    assert(!empty());
    return current;
  }

  void popFront() {
    assert(!empty());
    immutable block_size = stream.read!int();
    immutable refid = stream.read!int();
    immutable pos = stream.read!int();

    // void *p = pureMalloc(block_size-2*int.sizeof); // test for GC effectiveness
    data = new ubyte[block_size-2*int.sizeof];
    current = ReadBlob(refid,pos,stream.read(data));
    assert(current._data.ptr == data.ptr);
  }

  /// Returns a read if available. Otherwise null
  Nullable!ReadBlob read() {
    if (empty()) return Nullable!ReadBlob();
    auto prev = current;
    popFront();
    return current;
  }

  /// Returns the next matching read. Otherwise null
  ///
  /// Example:
  ///
  ///    auto current = ProcessReadBlob(stream.read_if!ProcessReadBlob((r) => !remove && r.is_mapped));
  Nullable!ReadBlob read_if(R)(bool delegate(R r) is_match) {
    while(!empty()) {
      read();
      if (is_match(R(current)))
        return current;
    }
    return Nullable!ReadBlob();
  }

}
