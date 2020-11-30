/*
    New style BAM reader. This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/

// This is a complete rewrite of Artem Tarasov's original reader.

module bio.std.experimental.hts.bam.reader;

import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.exception;
import std.file;
import std.format : format;
import std.stdio;
import std.string;
import std.typecons;
import std.bitmanip;

//TODO remove these dependecies 
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.constants;

import bio.core.utils.exception;

import bio.std.experimental.hts.bgzf;
import bio.std.experimental.hts.constants;

import bio.std.experimental.hts.bam.header;

template ReadFlags(alias flag) {
  // 0x01: template having multiple segments in sequencing
  @property bool is_paired()                nothrow { return cast(bool)(flag & 0x1); }
  // 0x2: Each segment properly aligned according to the aligner
  @property bool is_proper_pair()           nothrow { return cast(bool)(flag & 0x2); }
  // 0x4: Segment unmapped
  @property bool is_unmapped_raw()          nothrow { return cast(bool)(flag & 0x4); }
  @property bool is_mapped_raw()            nothrow { return cast(bool)(!(flag & 0x4)); }
  // 0x8: Next segment in template unmapped
  @property bool mate_is_unmapped()         nothrow { return cast(bool)(flag & 0x8); }
  // 0x10: SEQ being reverse complimented
  @property bool is_reverse_strand()        nothrow { return cast(bool)(flag & 0x10); }
  // 0x20: SEQ of the next segment in the template being reverse complemented
  @property bool mate_is_reverse_strand()   nothrow { return cast(bool)(flag & 0x20); }
  // 0x40: The first segment in the template
  @property bool is_first_of_pair()         nothrow { return cast(bool)(flag & 0x40); }
  // 0x80: The last segment in the template
  @property bool is_second_of_pair()        nothrow { return cast(bool)(flag & 0x80); }
  // 0x100: Secondary segment
  @property bool is_secondary_alignment()   nothrow { return cast(bool)(flag & 0x100); }
  // 0x200: Not passing filters, such as platform/vendor quality controls
  @property bool is_qc_fail() {
    assert(is_mapped_raw,to!string(this));
    return cast(bool)(flag & 0x200); }
  alias is_qc_fail failed_quality_control;
  /// 0x400: PCR or optical duplicate
  @property bool is_duplicate()             nothrow { return cast(bool)(flag & 0x400); }
  /// 0x800: Supplementary alignment
  @property bool is_supplementary()         nothrow { return cast(bool)(flag & 0x800); }
  @property string show_flags() {
    string res = format("b%b-%d",flag,flag);
    if (is_paired) res ~= " pair";
    if (is_proper_pair) res ~= " proper";
    if (is_mapped_raw) res ~= " mapped";
    if (is_unmapped_raw) res ~= " unmapped";
    if (mate_is_unmapped) res ~= " mate_unmapped";
    if (is_reverse_strand) res ~= " rev_strand";
    if (mate_is_reverse_strand) res ~= " mate_rev_strand";
    if (is_first_of_pair) res ~= " first_of_pair";
    if (is_second_of_pair) res ~= " second_of_pair";
    if (is_secondary_alignment) res ~= " secondary_aln";
    if (is_mapped_raw && is_qc_fail) res ~= " qc_fail";
    if (is_duplicate) res ~= " duplicate";
    if (is_supplementary) res ~= " suppl";
    return res;
  }
}

template CheckMapped(alias refid) {
  @property nothrow bool is_unmapped() {
    return is_unmapped_raw;
  }
  @property bool is_mapped() {
    debug {
      if (is_mapped_raw) {
        assert(refid != -1, "ref_id can not be -1 for mapped read");  // BAM spec
      }
    }
    return !is_unmapped_raw;
  }
}

enum Offset {
  bin_mq_nl=0, flag_nc=4, flag=6, l_seq=8, next_refID=12, next_pos=16, tlen=20, read_name=24
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

  mixin ReadFlags!(_flag);
  mixin CheckMapped!(refid);

  /*
  this(RefId ref_id, GenomePos read_pos, ubyte[] buf) {
    refid = ref_id;
    pos = read_pos;
    _data = buf;
  }
  */

  @property void cleanup() {
    destroy(_data);
    _data = null;
  }

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
  ushort _flag()          { return fetch!ushort(Offset.flag); }
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
    return "<** " ~ ReadBlob.stringof ~ " (data size " ~ to!string(_data.length) ~ ") " ~ to!string(refid) ~ ":" ~ to!string(pos) ~ " length " ~ to!string(sequence_length) ~ " flags " ~ show_flags() ~ ">";
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

  @property void cleanup() {
    _read2.cleanup;
  }

  @property nothrow bool isNull() {
    return _read2.isNull;
  }

  @property RefId ref_id() {
    enforce(_read2.is_mapped,"Trying to get ref_id an unmapped read " ~ to!string(_read2));
    return _read2.refid;
  }

  @property RefId raw_ref_id() {
    return _read2.refid;
  }

  @property nothrow uint _flag_nc() {
    return _read2._flag_nc;
  }

  @property nothrow ushort _flag() {
    return _read2._flag;
  }

  alias ref_id refid;

  /// Get the start position on the reference sequence (better use
  /// start_loc), i.e., the first base that gets consumed in the
  /// CIGAR.
  @property GenomePos start_pos() {
    assert(_read2.is_mapped,"Trying to get pos on an unmapped read"); // BAM spec
    asserte(_read2.pos < GenomePos.max);
    return cast(GenomePos)_read2.pos;
  }

  @property GenomePos raw_start_pos() {
    return cast(GenomePos)_read2.pos;
  }

  /// Get the end position on the reference sequence (better use end_loc)
  @property GenomePos end_pos() {
    assert(sequence_length > 0, "Trying to get end_pos on an empty read sequence");
    assert(!consumed_reference_bases.isNull);
    return start_pos + consumed_reference_bases;
  }

  @property GenomePos raw_end_pos() {
    return raw_start_pos + consumed_reference_bases;
  }

  @property GenomeLocation start_loc() {
    return GenomeLocation(ref_id,start_pos);
  }

  @property GenomeLocation end_loc() {
    return GenomeLocation(ref_id,end_pos);
  }

  @property @trusted MappingQuality mapping_quality() { // MAPQ
    assert(_read2.is_mapped,"Trying to get MAPQ on an unmapped read"); // BAM spec
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

  /// Count and caches consumed reference bases. Uses raw_cigar to
  /// avoid a heap allocation.
  @property @trusted Nullable!GenomePos consumed_reference_bases() {
    if (consumed_reference_bases2.isNull) {
      assert(_read2.is_mapped,"Trying to get consumed bases on an unmapped read"); // BAM spec
      assert(!read_name.isNull,"Trying to get CIGAR on RNAME is '*'"); // BAM spec
      auto raw = cast(uint[]) _read2.raw_cigar();
      if (raw.length==1 && raw[0] == '*')
        return consumed_reference_bases2; // null
      else {
        GenomePos bases = 0;
        for (size_t i = 0; i < raw.length; i++) {
          auto cigarop = CigarOperation(raw[i]);
          if (cigarop.is_reference_consuming)
            bases += cigarop.length;
        }
        consumed_reference_bases2 = bases;
      }
    }
    return consumed_reference_bases2;
  }

  /// Count query consumed bases. Uses raw_cigar to avoid a heap
  /// allocation.
  @property @trusted GenomePos consumed_query_bases() {
    assert(_read2.is_mapped,"Trying to get consumed bases on an unmapped read"); // BAM spec
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
      return bases;
    }
  }

  /// Return read name as a string. If unavailable returns
  /// null. Caches name.
  @property Nullable!string read_name() {
    if (read_name2.isNull) {
      assert(_read2.is_mapped,"Trying to get RNAME on an unmapped read"); // BAM spec
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
      assert(_read2.is_mapped,"Trying to get CIGAR on an unmapped read"); // BAM spec
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

  @property ubyte[] toBlob() {
    return _read2._data;
  }

  @property string posString() {
    return (is_mapped ? to!string(ref_id) ~ ":" ~ to!string(start_pos) : "unmapped");
  }

  string toString() {
    // return "<** " ~ ProcessReadBlob.stringof ~ ") " ~ to!string(_read2.refid) ~ ":" ~ to!string(_read2.pos) ~ " length " ~ to!string(sequence_length) ~ ">";
    return _read2.get.toString();
  }

}

/**
   BamReader2 is used for foreach loops
*/

struct BamReadBlobs {
  BgzfStream stream;
  BamHeader header;

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
  BamHeader header;
  Nullable!ReadBlob readbuf; // points to current read
  ubyte[] data; // in sync with readbuf

  this(string fn) {
    stream = BgzfStream(fn);
    fetch_bam_header(header, stream);
    if (!empty)
      popFront(); // preload front
  }

  bool empty() @property {
    return stream.eof();
  }

  // Returns first read available. If past eof returns null.
  Nullable!ReadBlob front() {
    return readbuf;
  }

  void popFront() {
    asserte(!empty); // should have been checked for
    immutable block_size = stream.read!int();
    immutable refid = stream.read!int();
    immutable pos = stream.read!int();

    // void *p = pureMalloc(block_size-2*int.sizeof); // test for GC effectiveness
    data = new ubyte[block_size-2*int.sizeof];
    readbuf = ReadBlob(refid,pos,stream.read(data));
    assert(readbuf._data.ptr == data.ptr);
  }

}

/**
   Reader - use on single thread only

   This one provides peek support. Peek looks one read ahead in the read stream.
*/

struct BamBlobReader {
  BgzfStream stream;
  BamHeader header;
  Nullable!ReadBlob peekbuf; // points to current read
  // ubyte[] data; // in sync with peekbuf

  this(string fn) {
    stream = BgzfStream(fn);
    fetch_bam_header(header, stream);
  }

  bool empty() @property {
    return peekbuf.isNull && stream.eof();
  }

  Nullable!ReadBlob peek() {
    if (peekbuf.isNull && !empty)
      fetch();
    return peekbuf;
  }

  /// Fetches the next read. If the peekbuf is not empty return that
  /// first and reset peekbuf.
  Nullable!ReadBlob fetch() {
    if (!peekbuf.isNull) {
      auto readbuf = peekbuf;
      peekbuf = Nullable!ReadBlob();
      return readbuf;
    }
    asserte(!empty); // should have been checked for
    immutable block_size = stream.read!int();
    immutable refid = stream.read!int();
    immutable pos = stream.read!int();

    // void *p = pureMalloc(block_size-2*int.sizeof); // test for GC effectiveness
    auto data = new ubyte[block_size-2*int.sizeof];
    peekbuf = ReadBlob(refid,pos,stream.read(data));
    return peekbuf;
  }

  /// Returns the next matching read. Otherwise null
  ///
  /// Example:
  ///
  ///    auto readbuf = ProcessReadBlob(stream.read_if!ProcessReadBlob((r) => !remove && r.is_mapped));
  /*
  Nullable!ReadBlob read_if(R)(bool delegate(R r) is_match) {
    while(!empty()) {
      auto readbuf = read();
      if (is_match(R(readbuf)))
        return readbuf;
      else
        return Nullable!ReadBlob();
    }
    return Nullable!ReadBlob();
  }
  */
}
