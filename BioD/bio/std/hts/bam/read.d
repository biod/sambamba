/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
/// $(P $(D BamRead) type provides convenient interface for working with SAM/BAM records.)
///
/// $(P All flags, tags, and fields can be accessed and modified.)
///
/// Examples:
/// ---------------------------
/// import std.conv;
/// ...
/// assert(!read.is_unmapped);              // check flag
/// assert(read.ref_id != -1);              // access field
///
/// int edit_distance = to!int(read["NM"]); // access tag
/// read["NM"] = 0;                         // modify tag
/// read["NM"] = null;                      // remove tag
/// read["NM"] = null;                      // no-op
///
/// foreach (tag, value; read)              // iterate over tags
///     writeln(tag, " ", value);           // and print their keys and values
///
/// read.sequence = "AGCAGACTACGTGTGCATAC"; // sets base qualities to 255
/// assert(read.base_qualities[0] == 255);
/// read.is_unmapped = true;                // set flag
/// read.ref_id = -1;                       // set field
/// ---------------------------
module bio.std.hts.bam.read;

import bio.core.base;
import bio.core.utils.format;

import bio.std.hts.bam.abstractreader;
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.writer;
import bio.std.hts.bam.tagvalue;
import bio.std.hts.bam.bai.bin;

import bio.std.hts.bam.md.core;

import bio.std.hts.utils.array;
import bio.std.hts.utils.value;
import bio.core.utils.switchendianness;

import bio.std.hts.thirdparty.msgpack : Packer, unpack;

import std.algorithm;
import std.range;
import std.conv;
import std.format;
import std.exception;
import std.system;
import std.traits;
import std.array;
import std.bitmanip;
import core.stdc.stdlib;

/**
  BAM record representation.
*/
struct BamRead {

    mixin TagStorage;

    /// Reference index in BAM file header
    @property    int ref_id()           const nothrow { return _refID; }
    /// ditto
    @property   void ref_id(int n)                    { _dup(); _refID = n; }

    /// Reference sequence name ('*' for unmapped reads)
    @property string ref_name()         const nothrow { return _ref_id_to_string(ref_id); }

    /// 0-based leftmost coordinate of the first matching base
    @property    int position()         const nothrow { return _pos; }
    /// ditto
    @property   void position(int n)                  { _dup(); _pos = n; _recalculate_bin(); }

    /// Indexing bin which this read belongs to. Recalculated when position is changed.
    @property    bio.std.hts.bam.bai.bin.Bin bin()              const nothrow { return Bin(_bin); }

    /// Mapping quality. Equals to 255 if not available, otherwise
    /// equals to rounded -10 * log10(P {mapping position is wrong}).
    @property  ubyte mapping_quality()  const nothrow { return _mapq; }
    /// ditto
    @property   void mapping_quality(ubyte n)         { _dup(); _mapq = n; }

    /// Flag bits (should be used on very rare occasions, see flag getters/setters below)
    @property ushort flag()             const nothrow { return _flag; }
    /// ditto
    @property   void flag(ushort n)                   { _dup(); _flag = n; }

    /// Sequence length. In fact, sequence.length can be used instead, but that might be
    /// slower if the compiler is not smart enough to optimize away unrelated stuff.
    @property    int sequence_length()  const nothrow { return _l_seq; }

    /// Mate reference ID
    @property    int mate_ref_id()      const nothrow { return _next_refID; }
    /// ditto
    @property   void mate_ref_id(int n)               { _dup(); _next_refID = n; }

    /// Mate reference sequence name ('*' for unmapped mates)
    @property string mate_ref_name()    const nothrow { return _ref_id_to_string(_next_refID); }

    /// Mate position
    @property    int mate_position()    const nothrow { return _next_pos; }
    /// ditto
    @property   void mate_position(int n)             { _dup(); _next_pos = n; }

    /// Template length
    @property    int template_length()  const nothrow { return _tlen; }
    /// ditto
    @property   void template_length(int n)           { _dup(); _tlen = n; }

    // ------------------------ FLAG GETTERS/SETTERS -------------------------------------- //

    /// Template having multiple segments in sequencing
    @property bool is_paired()                const nothrow { return cast(bool)(flag & 0x1); }
    /// ditto
    @property void is_paired(bool b)                { _setFlag( 0, b); }

    /// Each segment properly aligned according to the aligner
    @property bool proper_pair()              const nothrow { return cast(bool)(flag & 0x2); }
    /// ditto
    @property void proper_pair(bool b)              { _setFlag( 1, b); }

    /// Segment unmapped
    @property bool is_unmapped()              const nothrow { return cast(bool)(flag & 0x4); }
    /// ditto
    @property void is_unmapped(bool b)              { _setFlag( 2, b); }

    /// Next segment in the template unmapped
    @property bool mate_is_unmapped()         const nothrow { return cast(bool)(flag & 0x8); }
    /// ditto
    @property void mate_is_unmapped(bool b)         { _setFlag( 3, b); }

    /// Sequence being reverse complemented
    @property bool is_reverse_strand()        const nothrow { return cast(bool)(flag & 0x10); }
    /// ditto
    @property void is_reverse_strand(bool b)        { _setFlag( 4, b); }

    /// Sequence of the next segment in the template being reversed
    @property bool mate_is_reverse_strand()   const nothrow { return cast(bool)(flag & 0x20); }
    /// ditto
    @property void mate_is_reverse_strand(bool b)   { _setFlag( 5, b); }

    /// The first segment in the template
    @property bool is_first_of_pair()         const nothrow { return cast(bool)(flag & 0x40); }
    /// ditto
    @property void is_first_of_pair(bool b)         { _setFlag( 6, b); }

    /// The last segment in the template
    @property bool is_second_of_pair()        const nothrow { return cast(bool)(flag & 0x80); }
    /// ditto
    @property void is_second_of_pair(bool b)        { _setFlag( 7, b); }

    /// Secondary alignment
    @property bool is_secondary_alignment()   const nothrow { return cast(bool)(flag & 0x100); }
    /// ditto
    @property void is_secondary_alignment(bool b)   { _setFlag( 8, b); }

    /// Not passing quality controls
    @property bool failed_quality_control()   const nothrow { return cast(bool)(flag & 0x200); }
    /// ditto
    @property void failed_quality_control(bool b)   { _setFlag( 9, b); }

    /// PCR or optical duplicate
    @property bool is_duplicate()             const nothrow { return cast(bool)(flag & 0x400); }
    /// ditto
    @property void is_duplicate(bool b)             { _setFlag(10, b); }

    /// Supplementary alignment
    @property bool is_supplementary()         const nothrow { return cast(bool)(flag & 0x800); }
    /// ditto
    @property void is_supplementary(bool b)         { _setFlag(11, b); }

    /// Convenience function, returns '+' or '-' indicating the strand.
    @property char strand() const nothrow {
        return is_reverse_strand ? '-' : '+';
    }

    /// ditto
    @property void strand(char c) {
        enforce(c == '-' || c == '+', "Strand must be '-' or '+'");
        is_reverse_strand = c == '-';
    }

    /// Read name, length must be in 1..255 interval.
    @property string name() const nothrow {
        // notice -1: the string is zero-terminated, so we should strip that '\0'
        return cast(string)(_chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1]);
    }

    /// ditto
    @property void name(string new_name) {
        enforce(new_name.length >= 1 && new_name.length <= 255,
                "name length must be in 1-255 range");
        _dup();
        bio.std.hts.utils.array.replaceSlice(_chunk,
                 _chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1],
                 cast(ubyte[])new_name);
        _l_read_name = cast(ubyte)(new_name.length + 1);
    }

    /// List of CIGAR operations
    @property const(CigarOperation)[] cigar() const nothrow {
        return cast(const(CigarOperation)[])(_chunk[_cigar_offset .. _cigar_offset +
                                             _n_cigar_op * CigarOperation.sizeof]);
    }

    /// ditto
    @property void cigar(const(CigarOperation)[] c) {
        enforce(c.length < 65536, "Too many CIGAR operations, must be <= 65535");
        _dup();
        bio.std.hts.utils.array.replaceSlice(_chunk,
             _chunk[_cigar_offset .. _cigar_offset + _n_cigar_op * CigarOperation.sizeof],
             cast(ubyte[])c);

        _n_cigar_op = cast(ushort)(c.length);

        _recalculate_bin();
    }

    /// Extended CIGAR where M operators are replaced with =/X based
    /// on information from MD tag. Throws if the read doesn't have MD
    /// tag.
    auto extended_cigar() @property const {
        Value md = this["MD"];
        enforce(md.is_string);
        return makeExtendedCigar(cigar, mdOperations(*cast(string*)(&md)));
    }

    /// The number of reference bases covered by this read.
    /// $(BR)
    /// Returns 0 if the read is unmapped.
    int basesCovered() const {

        if (this.is_unmapped) {
            return 0; // actually, valid alignments should have empty cigar string
        }

        return reduce!"a + b.length"(0, filter!"a.is_reference_consuming"(cigar));
    }

    /// Human-readable representation of CIGAR string (same as in SAM format)
    string cigarString() const {
        char[] str;

        // guess size of resulting string
        str.reserve(_n_cigar_op * 3);

        foreach (cigar_op; cigar) {
            str ~= to!string(cigar_op.length);
            str ~= cigar_op.type;
        }
        return cast(string)str;
    }

    private @property inout(ubyte)[] raw_sequence_data() inout nothrow {
        return _chunk[_seq_offset .. _seq_offset + (_l_seq + 1) / 2];
    }

    /// Read-only random-access range for access to sequence data.
    static struct SequenceResult {

        private size_t _index;
        private ubyte[] _data = void;
        private size_t _len = void;
        private bool _use_first_4_bits = void;

        this(const(ubyte[]) data, size_t len, bool use_first_4_bits=true) {
            _data = cast(ubyte[])data;
            _len = len;
            _use_first_4_bits = use_first_4_bits;
        }

        ///
        @property bool empty() const {
            return _index >= _len;
        }

        ///
        @property bio.core.base.Base front() const {
            return opIndex(0);
        }

        ///
        @property bio.core.base.Base back() const {
            return opIndex(_len - 1);
        }

        /*
        I have no fucking idea why this tiny piece of code
        does NOT get inlined by stupid DMD compiler.

        Therefore I use string mixin instead.
        (hell yeah! Back to the 90s! C macros rulez!)

        private size_t _getActualPosition(size_t index) const
        {
            if (_use_first_4_bits) {
                // [0 1] [2 3] [4 5] [6 7] ...
                //            |
                //            V
                //   0     1     2     3
                return index >> 1;
            } else {
                // [. 0] [1 2] [3 4] [5 6] ...
                //            |
                //            V
                //   0     1     2     3
                return (index >> 1) + (index & 1);
            }
        }*/

        private static string _getActualPosition(string index) {
            return "((" ~ index ~") >> 1) + " ~
                   "(_use_first_4_bits ? 0 : ((" ~ index ~ ") & 1))";
        }

        private bool _useFirst4Bits(size_t index) const
        {
            auto res = index % 2 == 0;
            if (!_use_first_4_bits) {
                res = !res;
            }
            return res;
        }

        ///
        @property SequenceResult save() const {
            return SequenceResult(_data[mixin(_getActualPosition("_index")) .. $],
                                  _len - _index,
                                  _useFirst4Bits(_index));
        }

        ///
        SequenceResult opSlice(size_t i, size_t j) const {
            return SequenceResult(_data[mixin(_getActualPosition("_index + i")) .. $],
                                  j - i,
                                  _useFirst4Bits(_index + i));
        }

        ///
        @property bio.core.base.Base opIndex(size_t i) const {
            auto pos = _index + i;

            if (_use_first_4_bits)
            {
                if (pos & 1)
                    return Base.fromInternalCode(_data[pos >> 1] & 0xF);
                else
                    return Base.fromInternalCode(_data[pos >> 1] >> 4);
            }
            else
            {
                if (pos & 1)
                    return Base.fromInternalCode(_data[(pos >> 1) + 1] >> 4);
                else
                    return Base.fromInternalCode(_data[pos >> 1] & 0xF);
            }

            assert(false);
        }

        /// ditto
        @property void opIndexAssign(bio.core.base.Base base, size_t i) {
            auto pos = _index + i;

            if (_use_first_4_bits)
            {
                if (pos & 1)
                    _data[pos >> 1] &= 0xF0, _data[pos >> 1] |= base.internal_code;
                else
                    _data[pos >> 1] &= 0x0F, _data[pos >> 1] |= base.internal_code << 4;
            }
            else
            {
                if (pos & 1)
                    _data[(pos >> 1) + 1] &= 0x0F, _data[(pos >> 1) + 1] |= base.internal_code << 4;
                else
                    _data[pos >> 1] &= 0xF0, _data[pos >> 1] |= base.internal_code;
            }
        }


        ///
        void popFront() {
            ++_index;
        }

        ///
        void popBack() {
            --_len;
        }

        ///
        @property size_t length() const {
            return _len - _index;
        }

        alias length opDollar;

        void toString(scope void delegate(const(char)[]) dg) const {
            char[256] buf = void;
            size_t total = this.length;
            size_t written = 0;
            while (written < total) {
                size_t n = min(buf.length, total - written);
                foreach (j; 0 .. n)
                    buf[j] = opIndex(written + j).asCharacter;
                dg(buf[0 .. n]);
                written += n;
            }
        }
    }

    /// Random-access range of characters
    @property SequenceResult sequence() const {
        return SequenceResult(raw_sequence_data, sequence_length);
    }

    static assert(isRandomAccessRange!(ReturnType!sequence));

    /// Sets query sequence. Sets all base qualities to 255 (i.e. unknown).
    @property void sequence(string seq)
    {
        _dup();

        auto raw_length = (seq.length + 1) / 2;
        // set sequence
        auto replacement = uninitializedArray!(ubyte[])(raw_length + seq.length);
        replacement[raw_length .. $] = 0xFF;
        for (size_t i = 0; i < raw_length; ++i) {
            replacement[i] = cast(ubyte)(Base(seq[2 * i]).internal_code << 4);

            if (seq.length > 2 * i + 1)
                replacement[i] |= cast(ubyte)(Base(seq[2 * i + 1]).internal_code);
        }

        bio.std.hts.utils.array.replaceSlice(_chunk,
                     _chunk[_seq_offset .. _tags_offset],
                     replacement);

        _l_seq = cast(int)seq.length;
    }

    /// Quality data (phred-based scores)
    @property inout(ubyte)[] base_qualities() inout nothrow {
        return _chunk[_qual_offset .. _qual_offset + _l_seq * char.sizeof];
    }

    /// Set quality data - array length must be of the same length as the sequence.
    @property void base_qualities(const(ubyte)[] quality) {
        enforce(quality.length == _l_seq, "Quality data must be of the same length as sequence");
        _dup();
        _chunk[_qual_offset .. _qual_offset + _l_seq] = quality;
    }

    /*
      Constructs the struct from memory chunk
      */
    this(ubyte[] chunk, bool fix_byte_order=true) {

        // Switching endianness lazily is not a good idea:
        //
        // 1) switching byte order is pretty fast
        // 2) lazy switching for arrays can kill the performance,
        //    it has to be done once
        // 3) the code will be too complicated, whereas there're
        //    not so many users of big-endian systems
        //
        // In summa, BAM is little-endian format, so big-endian
        // users will suffer anyway, it's unavoidable.

        _chunk = chunk;
        this._is_slice = true;

        if (fix_byte_order && std.system.endian != Endian.littleEndian) {
            switchChunkEndianness();

            // Dealing with tags is the responsibility of TagStorage.
            fixTagStorageByteOrder();
        }
    }

    // Doesn't touch tags, only fields.
    // @@@TODO: NEEDS TESTING@@@
    private void switchChunkEndianness() {
        // First 8 fields are 32-bit integers:
        //
        // 0) refID                int
        // 1) pos                  int
        // 2) bin_mq_nl           uint
        // 3) flag_nc             uint
        // 4) l_seq                int
        // 5) next_refID           int
        // 6) next_pos             int
        // 7) tlen                 int
        // ----------------------------------------------------
        // (after them name follows which is string)
        //
        switchEndianness(_chunk.ptr, 8 * uint.sizeof);

        // Then we need to switch endianness of CIGAR data:
        switchEndianness(_chunk.ptr + _cigar_offset,
                         _n_cigar_op * uint.sizeof);
    }

    private size_t calculateChunkSize(string read_name,
                                      string sequence,
                                      in CigarOperation[] cigar)
    {
        return 8 * int.sizeof
                 + (read_name.length + 1) // tailing '\0'
                 + uint.sizeof * cigar.length
                 + ubyte.sizeof * ((sequence.length + 1) / 2)
                 + ubyte.sizeof * sequence.length;
    }

    /// Construct alignment from basic information about it.
    ///
    /// Other fields can be set afterwards.
    this(string read_name,                          // info for developers:
         string sequence,                           // these 3 fields are needed
         in CigarOperation[] cigar)                 // to calculate size of _chunk
    {
        enforce(read_name.length < 256, "Too long read name, length must be <= 255");
        enforce(cigar.length < 65536, "Too many CIGAR operations, must be <= 65535");

        if (this._chunk is null) {
            this._chunk = new ubyte[calculateChunkSize(read_name, sequence, cigar)];
        }

        this._refID      =  -1;         // set default values
        this._pos        =  -1;         // according to SAM/BAM
        this._mapq       = 255;         // specification
        this._next_refID =  -1;
        this._next_pos   =  -1;
        this._tlen       =   0;

        this._l_read_name = cast(ubyte)(read_name.length + 1); // tailing '\0'
        this._n_cigar_op  = cast(ushort)(cigar.length);
        this._l_seq       = cast(int)(sequence.length);

        // now all offsets can be calculated through corresponding properties

        // set default quality
        _chunk[_qual_offset .. _qual_offset + sequence.length] = 0xFF;

        // set CIGAR data
        auto _len = cigar.length * CigarOperation.sizeof;
        _chunk[_cigar_offset .. _cigar_offset + _len] = cast(ubyte[])(cigar);

        // set read_name
        auto _offset = _read_name_offset;
        _chunk[_offset .. _offset + read_name.length] = cast(ubyte[])read_name;
        _chunk[_offset + read_name.length] = cast(ubyte)'\0';

        this._is_slice = false;

        this.sequence = sequence;
    }

    /// Deep copy of the record.
    BamRead dup() @property const {
        BamRead result;
        result._chunk = this._chunk.dup;
        result._is_slice = false;
        result._modify_in_place = false;
        result._reader = cast()_reader;
        return result;
    }

    /// Compare two alignments, including tags
    /// (the tags must follow in the same order for equality).
    bool opEquals(BamRead other) const pure nothrow {
        // don't forget about _is_slice trick
        auto m = _cigar_offset;
        return _chunk[0 .. m - 1] == other._chunk[0 .. m - 1] &&
               _chunk[m .. $] == other._chunk[m .. $];
    }

    bool opEquals(const(BamRead) other) const pure nothrow {
        return opEquals(cast()other);
    }

    /// Size of the alignment record when output to stream in BAM format.
    /// Includes block_size as well (see SAM/BAM specification)
    @property size_t size_in_bytes() const {
        return int.sizeof + _chunk.length;
    }

    package void write(BamWriter writer) {
        writer.writeInteger(cast(int)(_chunk.length));

        ubyte old_byte = _chunk[_cigar_offset - 1];
        _chunk[_cigar_offset - 1] = 0;

        if (std.system.endian != Endian.littleEndian) {
            switchChunkEndianness();
            writer.writeByteArray(_chunk[0 .. _tags_offset]);
            switchChunkEndianness();
        } else {
            writer.writeByteArray(_chunk[0 .. _tags_offset]);
        }

        _chunk[_cigar_offset - 1] = old_byte;

        writeTags(writer);
    }

    /// Packs message in the following format:
    /// $(BR)
    /// MsgPack array with elements
    ///   $(OL
    ///     $(LI name - string)
    ///     $(LI flag - ushort)
    ///     $(LI reference sequence id - int)
    ///     $(LI leftmost mapping position (1-based) - int)
    ///     $(LI mapping quality - ubyte)
    ///     $(LI array of CIGAR operation lengths - int[])
    ///     $(LI array of CIGAR operation types - ubyte[])
    ///     $(LI mate reference sequence id - int)
    ///     $(LI mate position (1-based) - int)
    ///     $(LI template length - int)
    ///     $(LI segment sequence - string)
    ///     $(LI phred-base quality - ubyte[])
    ///     $(LI tags - map: string -> value))
    void toMsgpack(Packer)(ref Packer packer) const {
        packer.beginArray(13);
        packer.pack(cast(ubyte[])name);
        packer.pack(flag);
        packer.pack(ref_id);
        packer.pack(position + 1);
        packer.pack(mapping_quality);
        packer.pack(array(map!"a.length"(cigar)));
        packer.pack(array(map!"a.type"(cigar)));
        packer.pack(mate_ref_id);
        packer.pack(mate_position);
        packer.pack(template_length);
        packer.pack(to!string(sequence));
        packer.pack(base_qualities);

        packer.beginMap(tagCount());
        foreach (key, value; this) {
            packer.pack(key);
            packer.pack(value);
        }
    }

    /// String representation.
    /// $(BR)
    /// Possible formats are SAM ("%s") and JSON ("%j")
    void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const {
        if (size_in_bytes < 10000 && fmt.spec == 's') {
            auto p = cast(char*)alloca(size_in_bytes * 5);
            char* end = p;
            toSam(end);
            sink(p[0 .. end - p]);
        } else if (size_in_bytes < 5000 && fmt.spec == 'j') {
            auto p = cast(char*)alloca(size_in_bytes * 10 + 1000);
            char* end = p;
            toJson(end);
            sink(p[0 .. end - p]);
        } else if (fmt.spec == 's') {
            toSam(sink);
        } else if (fmt.spec == 'j') {
            toJson(sink);
        } else {
            throw new FormatException("unknown format specifier");
        }
    }

    /// ditto
    void toSam(Sink)(auto ref Sink sink) const
        if (isSomeSink!Sink)
    {
        sink.write(name);
        sink.write('\t');
        sink.write(flag);
        sink.write('\t');
        if (ref_id == -1 || _reader is null)
            sink.write('*');
        else
            sink.write(_reader.reference_sequences[ref_id].name);

        sink.write('\t');
        sink.write(position + 1);
        sink.write('\t');
        sink.write(mapping_quality);
        sink.write('\t');

        if (cigar.length == 0)
            sink.write('*');
        else
            foreach (op; cigar)
                op.toSam(sink);

        sink.write('\t');

        if (mate_ref_id == ref_id) {
            if (mate_ref_id == -1)
                sink.write("*\t");
            else
                sink.write("=\t");
        } else {
            if (mate_ref_id == -1 || _reader is null) {
                sink.write("*\t");
            } else {
                auto mate_name = _reader.reference_sequences[mate_ref_id].name;
                sink.write(mate_name);
                sink.write("\t");
            }
        }

        sink.write(mate_position + 1);
        sink.write('\t');
        sink.write(template_length);
        sink.write('\t');

        if (sequence_length == 0)
            sink.write('*');
        else
            foreach (char c; sequence)
                sink.write(c);
        sink.write('\t');

        if (base_qualities.length == 0 || base_qualities[0] == 0xFF)
            sink.write('*');
        else
            foreach (qual; base_qualities)
                sink.write(cast(char)(qual + 33));

        foreach (k, v; this) {
            sink.write('\t');
            sink.write(k);
            sink.write(':');
            v.toSam(sink);
        }
    }

    /// ditto
    string toSam()() const {
        return to!string(this);
    }

    /// JSON representation
    void toJson(Sink)(auto ref Sink sink) const
        if (isSomeSink!Sink)
    {
        sink.write(`{"qname":`); sink.writeJson(name);
        sink.write(`,"flag":`); sink.write(flag);

        sink.write(`,"rname":`);
        if (ref_id == -1 || _reader is null)
            sink.write(`"*"`);
        else
            sink.writeJson(_reader.reference_sequences[ref_id].name);

        sink.write(`,"pos":`); sink.write(position + 1);
        sink.write(`,"mapq":`); sink.write(mapping_quality);

        sink.write(`,"cigar":"`);
        if (cigar.empty)
            sink.write('*');
        else
            foreach (op; cigar)
                op.toSam(sink);
        sink.write('"');

        sink.write(`,"rnext":`);
        if (mate_ref_id == ref_id) {
            if (mate_ref_id == -1)
                sink.write(`"*"`);
            else
                sink.write(`"="`);
        } else if (mate_ref_id == -1 || _reader is null) {
            sink.write(`"*"`);
        } else {
            sink.writeJson(_reader.reference_sequences[mate_ref_id].name);
        }

        sink.write(`,"pnext":`); sink.write(mate_position + 1);
        sink.write(`,"tlen":`); sink.write(template_length);

        sink.write(`,"seq":"`);
        if (sequence_length == 0)
            sink.write('*');
        else
            foreach (char c; sequence)
                sink.write(c);
        sink.write('"');

        sink.write(`,"qual":`);
        sink.writeJson(base_qualities);

        sink.write(`,"tags":{`);

        bool not_first = false;
        foreach (k, v; this) {
            if (not_first)
                sink.write(',');
            sink.writeJson(k);
            sink.write(':');
            v.toJson(sink);
            not_first = true;
        }

        sink.write("}}");
    }

    /// ditto
    string toJson()() const {
        auto w = appender!(char[])();
        toJson((const(char)[] s) { w.put(s); });
        return cast(string)w.data;
    }

    /// Associates read with BAM reader. This is done automatically
    /// if this read is obtained through BamReader/Reference methods.
    void associateWithReader(bio.std.hts.bam.abstractreader.IBamSamReader reader) {
        _reader = reader;
    }

    /// Associated BAM/SAM reader.
    inout(bio.std.hts.bam.abstractreader.IBamSamReader) reader() @property inout {
        return _reader;
    }

    ///
    bool is_slice_backed() @property const {
        return _is_slice;
    }

    /// Raw representation of the read. Occasionally useful for dirty hacks!
    inout(ubyte)[] raw_data() @property inout {
        return _chunk;
    }

    /// ditto
    void raw_data(ubyte[] data) @property {
        _chunk = data;
    }

    package ubyte[] _chunk; // holds all the data,
                    // the access is organized via properties
                    // (see below)

private:

    // by specs, name ends with '\0'
    // let's use this byte for something useful!
    //
    // (Of course this places some restrictions on usage,
    //  but allows to reduce size of record.)
    bool _is_slice() @property const {
        return cast(bool)(_chunk[_cigar_offset - 1] & 1);
    }

    void _is_slice(bool is_slice) @property {
        _chunk[_cigar_offset - 1] &= 0b11111110;
        _chunk[_cigar_offset - 1] |= (is_slice ? 1 : 0);
    }

    // don't call _dup() if the record is modified
    bool _modify_in_place() @property const {
        return cast(bool)(_chunk[_cigar_offset - 1] & 2);
    }

    void _modify_in_place(bool in_place) @property {
        _chunk[_cigar_offset - 1] &= 0b11111101;
        _chunk[_cigar_offset - 1] |= (in_place ? 2 : 0);
    }

    IBamSamReader _reader;

    string _ref_id_to_string(int ref_id) const nothrow {
        if (_reader is null)
            return "?";
        if (ref_id < 0)
            return "*";
        return _reader.reference_sequences[ref_id].name;
    }

    // Official field names from SAM/BAM specification.
    // For internal use only
    @property  int _refID()      const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 0));
    }

    @property  int _pos()        const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 1));
    }

    @property uint _bin_mq_nl()  const nothrow pure @system {
        return *(cast(uint*)(_chunk.ptr + int.sizeof * 2));
    }

    @property uint _flag_nc()    const nothrow {
        return *(cast(uint*)(_chunk.ptr + int.sizeof * 3));
    }

    @property  int _l_seq()      const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 4));
    }

    @property  int _next_refID() const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 5));
    }

    @property  int _next_pos()   const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 6));
    }

    @property  int _tlen()       const nothrow {
        return *(cast( int*)(_chunk.ptr + int.sizeof * 7));
    }

    // Setters, also only for internal use
    @property void _refID(int n)       { *(cast( int*)(_chunk.ptr + int.sizeof * 0)) = n; }
    @property void _pos(int n)         { *(cast( int*)(_chunk.ptr + int.sizeof * 1)) = n; }
    @property void _bin_mq_nl(uint n)  { *(cast(uint*)(_chunk.ptr + int.sizeof * 2)) = n; }
    @property void _flag_nc(uint n)    { *(cast(uint*)(_chunk.ptr + int.sizeof * 3)) = n; }
    @property void _l_seq(int n)       { *(cast( int*)(_chunk.ptr + int.sizeof * 4)) = n; }
    @property void _next_refID(int n)  { *(cast( int*)(_chunk.ptr + int.sizeof * 5)) = n; }
    @property void _next_pos(int n)    { *(cast( int*)(_chunk.ptr + int.sizeof * 6)) = n; }
    @property void _tlen(int n)        { *(cast( int*)(_chunk.ptr + int.sizeof * 7)) = n; }

    // Additional useful properties, also from SAM/BAM specification
    //
    //             The layout of bin_mq_nl and flag_nc is as follows
    //                     (upper bits -------> lower bits):
    //
    // bin_mq_nl [ { bin (16b) }  { mapping quality (8b) } { read name length (8b) } ]
    //
    // flag_nc   [ { flag (16b) } { n_cigar_op (16b) } ]
    //
    @property ushort _bin()         const nothrow {
        return _bin_mq_nl >> 16;
    }
    @property  ubyte _mapq()        const nothrow {
        return (_bin_mq_nl >> 8) & 0xFF;
    }
    @property  ubyte _l_read_name() const nothrow pure {
        return _bin_mq_nl & 0xFF;
    }
    @property ushort _flag()        const nothrow {
        return _flag_nc >> 16;
    }
    @property ushort _n_cigar_op()  const nothrow {
        return _flag_nc & 0xFFFF;
    }

    // Setters for those properties
    @property void _bin(ushort n)         { _bin_mq_nl = (_bin_mq_nl &  0xFFFF) | (n << 16); }
    @property void _mapq(ubyte n)         { _bin_mq_nl = (_bin_mq_nl & ~0xFF00) | (n << 8); }
    @property void _l_read_name(ubyte n)  { _bin_mq_nl = (_bin_mq_nl & ~0xFF  ) | n; }
    @property void _flag(ushort n)        { _flag_nc   = (_flag_nc   &  0xFFFF) | (n << 16); }
    @property void _n_cigar_op(ushort n)  { _flag_nc   = (_flag_nc   & ~0xFFFF) | n; }

    // Offsets of various arrays in bytes.
    // Currently, are computed each time, so if speed will be an issue,
    // they can be made fields instead of properties.
    @property size_t _read_name_offset() const nothrow pure {
        return 8 * int.sizeof;
    }

    @property size_t _cigar_offset()     const nothrow pure {
        return _read_name_offset + _l_read_name * char.sizeof;
    }

    @property size_t _seq_offset()       const nothrow {
        return _cigar_offset + _n_cigar_op * uint.sizeof;
    }

    @property size_t _qual_offset()      const nothrow {
        return _seq_offset + (_l_seq + 1) / 2;
    }

    // Offset of auxiliary data
    @property size_t _tags_offset()      const nothrow {
        return _qual_offset + _l_seq;
    }

    // Sets n-th flag bit to boolean value b.
    void _setFlag(int n, bool b) {
        assert(n < 16);
        // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
        ushort mask = cast(ushort)(1 << n);
        _flag = (_flag & ~cast(int)(mask)) | ((-cast(int)b) & mask);
    }

    // If _chunk is still a slice, not an array, duplicate it.
    // Used when some part of alignment record is modified by user.
    //
    // Basically, it's sort of copy-on-write: a lot of read-only alignments
    // may point to the same location, but every modified one allocates its
    // own chunk of memory.
    void _dup() {
        if (_is_slice && !_modify_in_place) {
            _chunk = _chunk.dup;
            _is_slice = false;
        }
    }

public:
    // Calculates bin number.
    void _recalculate_bin() {
        _bin = reg2bin(position, position + basesCovered());
    }
}


/// Lazy tag storage.
///
///   Provides hash-like access and opportunity to iterate
///   storage like an associative array.
mixin template TagStorage() {

    // Provides access to chunk of memory which contains tags.
    // This way, every time _tags_offset gets updated
    // (due to update of cigar string/read name/sequence and memory move),
    // the change is reflected automatically in tag storage.
    private @property const(ubyte)[] _tags_chunk() const  {
        return _chunk[_tags_offset .. $];
    }

    /// Hash-like access to tags. Time complexity is $(BIGOH number of tags).
    /// $(BR)
    /// If tag with such $(I key) is not found, returned value 'is nothing'.
    /// $(BR)
    /// If key length is different from 2, exception is thrown.
    /// $(BR)
    /// Special case when $(I value) represents nothing is used for removing tag
    /// (assuming that no more than one with this key is presented in the record).
    ///
    /// Examples:
    /// ----------------------------
    /// auto v = read["NM"];
    /// assert(v.is_integer);
    ///
    /// auto v = read["MN"];
    /// assert(v.is_nothing); // no such tag
    ///
    /// read["NM"] = 3; // converted to bio.std.hts.bam.tagvalue.Value implicitly
    ///
    /// read["NM"] = null; // removes tag
    /// assert(read["NM"].is_nothing);
    /// ----------------------------
    bio.std.hts.bam.tagvalue.Value opIndex(string key) const {
        enforce(key.length == 2, "Key length must be 2");
        auto __tags_chunk = _tags_chunk; // _tags_chunk is evaluated lazily
        if (__tags_chunk.length < 4)
            return Value(null);

       size_t offset = 0;
       while (offset + 1 < __tags_chunk.length) {
           if (__tags_chunk[offset .. offset + 2] == key) {
               offset += 2;
               return readValue(offset, __tags_chunk);
           } else {
               offset += 2;
               skipValue(offset, __tags_chunk);
           }
       }
       return Value(null);
    }

    /// ditto
    void opIndexAssign(T)(T value, string key)
        if (is(T == Value) || __traits(compiles, GetTypeId!T))
    {
        static if(is(T == Value)) {
            enforce(key.length == 2, "Key length must be 2");
            auto __tags_chunk = _tags_chunk;

            _dup();

            size_t offset = 0;
            while (offset + 1 < __tags_chunk.length) {
                if (__tags_chunk[offset .. offset + 2] == key) {
                    if (value.is_nothing) {
                        // special case - remove tag
                        removeValueAt(offset);
                    } else {
                        replaceValueAt(offset + 2, value);
                    }
                    return;
                } else {
                    offset += 2;
                    skipValue(offset, __tags_chunk);
                }
            }

            if (!value.is_nothing)
                appendTag(key, value);
        } else {
            opIndexAssign(Value(value), key);
        }
    }

    /// Append new tag to the end, skipping check if it already exists. $(BIGOH 1)
    void appendTag(string key, Value value) {
        auto oldlen = _chunk.length;
        _chunk.length = _chunk.length + sizeInBytes(value) + 2 * char.sizeof;
        _chunk[oldlen .. oldlen + 2] = cast(ubyte[])key;
        emplaceValue(_chunk.ptr + oldlen + 2, value);
    }

    /// Remove all tags
    void clearAllTags() {
        _chunk.length = _tags_offset;
    }

    /// Number of tags. $(BIGOH number of tags)
    size_t tagCount() {
        size_t result = 0;
        size_t offset = 0;
        auto __tags_chunk = _tags_chunk;
        while (offset + 1 < __tags_chunk.length) {
            offset += 2;
            skipValue(offset, __tags_chunk);
            result += 1;
        }
        return result;
    }

    // replace existing tag
    private void replaceValueAt(size_t offset, Value value) {
        // offset points to the beginning of the value
        auto begin = offset;
        auto __tags_chunk = _tags_chunk;
        skipValue(offset, __tags_chunk); // now offset is updated and points to the end
        auto end = offset;

        prepareSlice(_chunk, __tags_chunk[begin .. end], sizeInBytes(value));

        emplaceValue(_chunk.ptr + _tags_offset + begin, value);
    }

    // remove existing tag
    private void removeValueAt(size_t begin) {
        // offset points to the beginning of the value
        auto offset = begin + 2;
        auto __tags_chunk = _tags_chunk;
        skipValue(offset, __tags_chunk);
        auto end = offset;
        // this does the job (see prepareSlice code)
        prepareSlice(_chunk, __tags_chunk[begin .. end], 0);
    }

    ///  Provides opportunity to iterate over tags.
    int opApply(scope int delegate(const ref string k, const ref Value v) dg) const {
        size_t offset = 0;
        auto __tags_chunk = _tags_chunk;
        while (offset + 1 < __tags_chunk.length) {
            auto key = cast(string)__tags_chunk[offset .. offset + 2];
            offset += 2;
            auto val = readValue(offset, __tags_chunk);
            auto res = dg(key, val);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

    /// Returns the number of tags. Time complexity is $(BIGOH number of tags)
    size_t tagCount() const {
        size_t res = 0;
        size_t offset = 0;
        auto __tags_chunk = _tags_chunk;
        while (offset + 1 < __tags_chunk.length) {
            offset += 2;
            skipValue(offset, __tags_chunk);
            res += 1;
        }
        return res;
    }

    private void writeTags(BamWriter writer) {
        if (std.system.endian == Endian.littleEndian) {
            writer.writeByteArray(_tags_chunk[]);
        } else {
            fixTagStorageByteOrder();
            writer.writeByteArray(_tags_chunk[]);
            fixTagStorageByteOrder();
        }
    }

    // Reads value which starts from (_tags_chunk.ptr + offset) address,
    // and updates offset to the end of value. O(1)
    private Value readValue(ref size_t offset, const(ubyte)[] tags_chunk) const {
        char type = cast(char)tags_chunk[offset++];
        return readValueFromArray(type, tags_chunk, offset);
    }

    // Increases offset so that it points to the next value. O(1).
    private void skipValue(ref size_t offset, const(ubyte)[] tags_chunk) const {
        char type = cast(char)tags_chunk[offset++];
        if (type == 'Z' || type == 'H') {
            while (tags_chunk[offset++] != 0) {}
        } else if (type == 'B') {
            char elem_type = cast(char)tags_chunk[offset++];
            auto length = *(cast(uint*)(tags_chunk.ptr + offset));
            offset += uint.sizeof + charToSizeof(elem_type) * length;
        } else {
            offset += charToSizeof(type);
        }
    }

    /*
      Intended to be used in constructor for initial endianness fixing
      in case the library is used on big-endian system.

      NOT TESTED AT ALL!!!
    */
    private void fixTagStorageByteOrder() {
        /* TODO: TEST ON BIG-ENDIAN SYSTEM!!! */
        const(ubyte)* p = _tags_chunk.ptr;
        const(ubyte)* end = p + _chunk.length;
        while (p < end) {
            p += 2; // skip tag name
            char type = *(cast(char*)p);
            ++p; // skip type
            if (type == 'Z' || type == 'H') {
                while (*p != 0) { // zero-terminated
                    ++p;          // string
                }
                ++p; // skip '\0'
            } else if (type == 'B') { // array
                char elem_type = *(cast(char*)p);
                uint size = charToSizeof(elem_type);
                switchEndianness(p, uint.sizeof);
                uint length = *(cast(uint*)p);
                p += uint.sizeof; // skip length
                if (size != 1) {
                    for (auto j = 0; j < length; j++) {
                        switchEndianness(p, size);
                        p += size;
                    }
                } else {
                    // skip
                    p += length;
                }
            } else {
                uint size = charToSizeof(type);
                if (size != 1) {
                    switchEndianness(p, size);
                    p += size;
                } else {
                    ++p;
                }
            }
        }
    }
}

unittest {

    import std.algorithm;
    import std.stdio;
    import std.math;

    // stderr.writeln("Testing BamRead behaviour...");
    auto read = BamRead("readname",
                        "AGCTGACTACGTAATAGCCCTA",
                        [CigarOperation(22, 'M')]);
    assert(read.sequence_length == 22);
    assert(read.cigar.length == 1);
    assert(read.cigarString() == "22M");
    assert(read.name == "readname");
    assert(equal(read.sequence(), "AGCTGACTACGTAATAGCCCTA"));

    read.name = "anothername";
    assert(read.name == "anothername");
    assert(read.cigarString() == "22M");

    read.base_qualities = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                           13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
    assert(reduce!"a+b"(0, read.base_qualities) == 253);

    read["RG"] = 15;
    assert(read["RG"] == 15);

    read["X1"] = [1, 2, 3, 4, 5];
    assert(read["X1"] == [1, 2, 3, 4, 5]);

    read.cigar = [CigarOperation(20, 'M'), CigarOperation(2, 'X')];
    assert(read.cigarString() == "20M2X");

    read["RG"] = cast(float)5.6;
    assert(approxEqual(to!float(read["RG"]), 5.6));

    read.sequence = "AGCTGGCTACGTAATAGCCCT";
    assert(read.sequence_length == 21);
    assert(read.base_qualities.length == 21);
    assert(read.base_qualities[20] == 255);
    assert(equal(read.sequence(), "AGCTGGCTACGTAATAGCCCT"));
    assert(retro(read.sequence)[2] == 'C');
    assert(retro(read.sequence)[0] == 'T');
    assert(read.sequence[4] == 'G');
    assert(read.sequence[0] == 'A');
    assert(equal(read.sequence[0..8], "AGCTGGCT"));
    assert(equal(read.sequence[3..5], "TG"));
    assert(equal(read.sequence[3..9][1..4], "GGC"));

    read["X1"] = 42;
    assert(read["X1"] == 42);

    assert(read.tagCount() == 2);

    read["X1"] = null;
    assert(read["X1"].is_nothing);
    assert(read.tagCount() == 1);
    read.sequence = "GTAAGCTGGCACTAGCAGCCT";
    read.cigar = [CigarOperation(read.sequence_length, 'M')];
    read["RG"] = null;
    read["RG"] = "readgroup1";
    assert(read.tagCount() == 1);
    read["RG"] = null;
    assert(read.tagCount() == 0);

    read.sequence[5] = Base('N');
    read.sequence[6] = Base('A');
    read.sequence[7] = Base('C');
    read.sequence[8] = Base('G');
    read.base_qualities[5] = 42;
    assert(read.sequence[5 .. 9].equal("NACG"));
    assert(read.base_qualities[5] == 42);

    // Test MsgPack serialization/deserialization

    {
    import std.typecons;
    static import bio.std.hts.thirdparty.msgpack;
    auto packer = bio.std.hts.thirdparty.msgpack.packer(Appender!(ubyte[])());
    read.toMsgpack(packer);
    auto data = packer.stream.data;
    auto rec = unpack(data).via.array;
    assert(rec[0].as!(ubyte[]) == "anothername");
    assert(rec[5].as!(int[]) == [21]);
    assert(rec[6].as!(ubyte[]) == ['M']);
    assert(rec[10].as!(ubyte[]) == to!string(read.sequence));
    }

    read.clearAllTags();
    assert(read.tagCount() == 0);
}

/// $(P BamRead wrapper which precomputes $(D end_position) = $(D position) + $(D basesCovered()).)
///
/// $(P Computation of basesCovered() takes quite a few cycles. Therefore in places where this
/// property is frequently accessed, it makes sense to precompute it for later use.)
///
/// $(P The idea is that this should be a drop-in replacement for BamRead in algorithms,
/// as the struct uses 'alias this' construction for the wrapped read.)
struct EagerBamRead(R=BamRead) {
    ///
    this(R read) {
        this.read = read;
        this.end_position = read.position + read.basesCovered();
    }

    ///
    R read;
    ///
    alias read this;

    /// End position on the reference, computed as position + basesCovered().
    int end_position;

    ///
    EagerBamRead dup() @property const {
        return EagerBamRead(read.dup);
    }
}

static assert(is(EagerBamRead!BamRead : BamRead));

/// Checks if $(D T) behaves like $(D BamRead)
template isBamRead(T)
{
    static if (is(Unqual!T : BamRead))
        enum isBamRead = true;
    else
        enum isBamRead = __traits(compiles,
        {
            T t; bool p;
            p = t.ref_id == 1;          p = t.position == 2;          p = t.bin.id == 3;
            p = t.mapping_quality == 4; p = t.flag == 5;              p = t.sequence_length == 6;
            p = t.mate_ref_id == 7;     p = t.mate_position == 8;     p = t.template_length == 9;
            p = t.is_paired;            p = t.proper_pair;            p = t.is_unmapped;
            p = t.mate_is_unmapped;     p = t.mate_is_reverse_strand; p = t.is_first_of_pair;
            p = t.is_second_of_pair;    p = t.is_secondary_alignment; p = t.failed_quality_control;
            p = t.is_duplicate;         p = t.strand == '+';          p = t.name == "";
            p = t.cigar[0].type == 'M'; p = t.basesCovered() > 42;    p = t.cigarString() == "";
            p = t.sequence[0] == 'A';   p = t.base_qualities[0] == 0;
        });
}

// Comparison function for 'queryname' sorting order based on https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/SAMRecordQueryNameComparator.java
bool compareReadNamesAsPicard(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isBamRead!R2)
 {
    if(a1.name == a2.name)
    {
        if(a1.is_paired() || a2.is_paired())
        {
            if(!a1.is_paired())
                return false;
            if(!a2.is_paired())
                return true;

            if(a1.is_first_of_pair() && a2.is_second_of_pair())
                return true;

            if(a1.is_second_of_pair() && a2.is_first_of_pair())
                return false;
                
        }

        if(a1.strand() != a2.strand())
        {
            return a1.strand() == '-' ? false : true;
        }

        if(a1.is_secondary_alignment() != a2.is_secondary_alignment())
        {
            return a2.is_secondary_alignment();
        }

        if(a1.is_supplementary() != a2.is_supplementary())
        {
            return a2.is_supplementary();
        }

        if(!a1["HI"].is_nothing)
        {
                if(a2["HI"].is_nothing)
                        return true;

                int i1 = to!int(a1["HI"]);
                int i2 = to!int(a2["HI"]);
                return i1 < i2;
        }
        else
        if(!a2["HI"].is_nothing)
                return false;
    }
    return a1.name < a2.name;
}

bool compareReadNamesAsPicard(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isSomeString!R2)
{
    return a1.name < a2;
}

bool compareReadNamesAsPicard(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isSomeString!R1 && isBamRead!R2)
{
    return a1 < a2.name;
}

/// $(P Comparison function for 'queryname' sorting order
/// (return whether first read is 'less' than second))
///
/// $(P This function can be called on:
///   $(UL
///     $(LI two reads)
///     $(LI read and string in any order)))
bool compareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isBamRead!R2)
{
    return a1.name < a2.name;
}

bool compareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isSomeString!R2)
{
    return a1.name < a2;
}

bool compareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isSomeString!R1 && isBamRead!R2)
{
    return a1 < a2.name;
}

int mixedStrCompare(string a, string b) {
  import std.ascii : isDigit;
  while (!a.empty && !b.empty) {
    if (a.front.isDigit && b.front.isDigit) {
      // skip zeros
      int za, zb;
      while (!a.empty && a.front == '0') { ++za; a.popFront(); }
      while (!b.empty && b.front == '0') { ++zb; b.popFront(); }

      // skip equal digits
      while (!a.empty && !b.empty && a.front.isDigit && a.front == b.front) {
        a.popFront();
        b.popFront();
      }

      if (!a.empty && !b.empty && a.front.isDigit && b.front.isDigit) {
        // the number of leading digits in each string is non-zero
        size_t i = 0, maxi = min(a.length, b.length);
        while (i < maxi && a[i].isDigit && b[i].isDigit) ++i;
        if (i < a.length && a[i].isDigit) return 1; // a contains more digits
        if (i < b.length && b[i].isDigit) return -1; // b contains more digits
        // the counts are equal, compare first digits
        return cast(byte)a.front - cast(byte)b.front;
      } else if (!a.empty && a.front.isDigit) return 1;
        else if (!b.empty && b.front.isDigit) return -1;
        else if (za != zb) return za - zb;  // order by the number of leading zeros
    } else {
      // lexicographical comparison for non-digits
      if (a.front != b.front) return cast(byte)a.front - cast(byte)b.front;
      a.popFront(); b.popFront();
    }
  }
  return (!a.empty) ? 1 : (!b.empty) ? -1 : 0;
}

/// $(P Comparison function for 'queryname' sorting order as in Samtools
/// (returns whether first read is 'less' than second in a 'mixed' order,
///  i.e. numbers inside the strings are compared by their integer value))
///
/// $(P This function can be called on:
///   $(UL
///     $(LI two reads)
///     $(LI read and string in any order)))
bool mixedCompareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isBamRead!R2)
{
    return mixedStrCompare(a1.name, a2.name) < 0;
}

bool mixedCompareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isSomeString!R2)
{
    return mixedStrCompare(a1.name, a2) < 0;
}

bool mixedCompareReadNames(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isSomeString!R1 && isBamRead!R2)
{
    return mixedStrCompare(a1, a2.name) < 0;
}

unittest {
  assert(mixedStrCompare("BC0123", "BC01234") < 0);
  assert(mixedStrCompare("BC0123", "BC0123Z") < 0);
  assert(mixedStrCompare("BC01234", "BC01234") == 0);
  assert(mixedStrCompare("BC0123DEF45", "BC01234DEF45") < 0);
  assert(mixedStrCompare("BC01236DEF45", "BC01234DEF45") > 0);
  assert(mixedStrCompare("BC012", "BC0012") < 0);
  assert(mixedStrCompare("BC0012DE0034", "BC0012DE34") > 0);
  assert(mixedStrCompare("BC12DE0034", "BC012DE34") < 0);
  assert(mixedStrCompare("1235", "1234") > 0);
}

// small utility function to get the value of the HI tag and return '0' 
// if it is not defined
private int getHI(R)(auto ref R r) 
    if (isBamRead!R)
{
        auto v = r["HI"];
        if (v.is_nothing)
            return 0;
        return to!int(v);
}

/// $(P Comparison function for 'queryname' sorting order setting mates
/// of the same alignments adjacent with the first mate coming before
/// the second mate)
bool compareReadNamesAndMates(R1, R2)(const auto ref R1 r1, const auto ref R2 r2) 
    if (isBamRead!R1 && isBamRead!R2)
{
    if (r1.name == r2.name) {
        if (getHI(r1) == getHI(r2))
            return r1.flag() < r2.flag();
        return getHI(r1) < getHI(r2);
    }
    return r1.name < r2.name;
}

/// $(P Comparison function for 'queryname' sorting order as in Samtools
/// setting mates of the same alignments adjacent with the first mate 
/// coming before the second mate)
bool mixedCompareReadNamesAndMates(R1, R2)(const auto ref R1 r1, const auto ref R2 r2)
    if (isBamRead!R1 && isBamRead!R2) 
{
    if (mixedStrCompare(r1.name, r2.name) == 0) {
        if (getHI(r1) == getHI(r2))
            return r1.flag() < r2.flag();
        return getHI(r1) < getHI(r2);    
    }
    return mixedStrCompare(r1.name, r2.name) < 0;
}

/// $(P Comparison function for 'coordinate' sorting order
/// (returns whether first read is 'less' than second))
///
/// $(P This function can be called on:
///   $(UL
///     $(LI two reads (in this case, reference IDs are also taken into account))
///     $(LI read and integer in any order)))

/// This function takes read direction into account (used for original samtools style sorting)
bool compareCoordinatesAndStrand(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isBamRead!R2)
{
    if (a1.ref_id == -1) return false; // unmapped reads should be last
    if (a2.ref_id == -1) return true;
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    if (a1.position < a2.position) return true;
    if (a1.position > a2.position) return false;
    return !a1.is_reverse_strand && a2.is_reverse_strand;
}

bool compareCoordinates(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isBamRead!R2)
{
    if (a1.ref_id == -1) return false; // unmapped reads should be last
    if (a2.ref_id == -1) return true;
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    return (a1.position < a2.position);
}

bool compareCoordinates(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isBamRead!R1 && isIntegral!R2)
{
    return a1.position < a2;
}

bool compareCoordinates(R1, R2)(const auto ref R1 a1, const auto ref R2 a2)
    if (isIntegral!R1 && isBamRead!R2)
{
    return a1 < a2.position;
}

static assert(isTwoWayCompatible!(compareReadNames, BamRead, string));
static assert(isTwoWayCompatible!(compareCoordinates, BamRead, int));

/// Allows modification of the read in-place even if it's slice-backed.
struct UniqueRead(R) {
    R read;
    alias read this;

    this(R read) {
        this.read = read;
        this.read._modify_in_place = true;
    }

    ~this() {
        this.read._modify_in_place = false;
    }
}

/// ditto
auto assumeUnique(R)(auto ref R read) if (isBamRead!R) {
    return UniqueRead!R(read);
}
