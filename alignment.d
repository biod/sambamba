module alignment;

import tagvalue;
private import bai.bin;

import std.algorithm;
import std.range;
import std.conv;
import std.exception;
import std.stream;
import std.system;

import utils.array;
import utils.value;
import utils.switchendianness;

/**
  Represents single CIGAR operation
 */
struct CigarOperation {
    static assert(CigarOperation.sizeof == uint.sizeof);
    /*
        WARNING!

      It is very essential that the size of 
      this struct is EXACTLY equal to uint.sizeof!

      The reason is to avoid copying of arrays during alignment parsing.

      Namely, when some_pointer points to raw cigar data,
      we can just do a cast. This allows to access those data
      directly, not doing any memory allocations. 
    */

    private uint raw; /// raw data from BAM

    private static ubyte char2op(char c) {
        switch(c) {
            case 'M': return 0;
            case 'I': return 1;
            case 'D': return 2;
            case 'N': return 3;
            case 'S': return 4;
            case 'H': return 5;
            case 'P': return 6;
            case '=': return 7;
            case 'X': return 8;
            default:  return 8; // FIXME
        }
    }

    this(uint length, char operation) {
        enforce(length < (1<<28), "Too big length of CIGAR operation");
        raw = (length << 4) | char2op(operation);
    }

    /// operation length
    uint length() @property const nothrow {
        return raw >> 4;
    }
  
    /// CIGAR operation as one of MIDNSHP=X
    char operation() @property const nothrow {
        if ((raw & 0xF) > 8) {
            /* FIXME: what to do in this case? */
            return '\0';
        }

        return "MIDNSHP=X"[raw & 0xF];
    }
}

/** 
  Represents single read
*/
struct Alignment {

    mixin TagStorage;

    /// This is an absolutely awful hack which allows you to write
    /// alignment.tags["X0"] / foreach(k, v; alignment.tags)
    /// instead of alignment["X0"] / foreach(k, v; alignment)
    Alignment tags() @property {
        return this;
    }

    // /////////////////////////////////////////////////////////////////////////
    //
    // Layout of Alignment in memory:
    //
    // TagStorage      <- tags
    //     ubyte[]     <- slice of the whole chunk. Must maintain it
    //         size_t                               because struct cannot hold
    //         ubyte*                               a pointer to itself in D.
    //     size_t      <- tags offset in the chunk
    // ubyte[]         <- _chunk: binary representation of alignment, as in BAM.
    //     size_t
    //     ubyte*
    // bool            <- _is_slice: indicates whether _chunk is a slice or not
    //
    // //////////////////////////////////////////////////////////////////////////

    @property    int ref_id()           const nothrow { return _refID; }
    @property    int position()         const nothrow { return _pos; }
    @property ushort bin()              const nothrow { return _bin; }
    @property  ubyte mapping_quality()  const nothrow { return _mapq; }
    @property ushort flag()             const nothrow { return _flag; }
    @property    int sequence_length()  const nothrow { return _l_seq; }
    @property    int next_ref_id()      const nothrow { return _next_refID; }
    @property    int next_pos()         const nothrow { return _next_pos; }
    @property    int template_length()  const nothrow { return _tlen; }

    /// Set reference id.
    @property void ref_id(int n)              { _dup(); _refID = n; }
    /// Sets 0-based leftmost coordinate. Bin is automatically recalculated.
    @property void position(int n)            { _dup(); _pos = n; _recalculate_bin(); }
    /// Set mapping quality
    @property void mapping_quality(ubyte n)   { _dup(); _mapq = n; }
    /// Set flag 
    @property void flag(ushort n)             { _dup(); _flag = n; }
    /// Set mate reference id.
    @property void next_ref_id(int n)         { _dup(); _next_refID = n; }
    /// Set mate position
    @property void next_pos(int n)            { _dup(); _next_pos = n; }
    /// Set template length
    @property void template_length(int n)     { _dup(); _tlen = n; }
  
    /// Template having multiple segments in sequencing
    @property bool is_paired()                const nothrow { return cast(bool)(flag & 0x1); }
    /// Each segment properly aligned according to the aligner
    @property bool proper_pair()              const nothrow { return cast(bool)(flag & 0x2); }
    /// Segment unmapped
    @property bool is_unmapped()              const nothrow { return cast(bool)(flag & 0x4); }
    /// Next segment in the template unmapped
    @property bool mate_is_unmapped()         const nothrow { return cast(bool)(flag & 0x8); }
    /// Sequence being reverse complemented
    @property bool is_reverse_strand()        const nothrow { return cast(bool)(flag & 0x10); }
    /// Sequence of the next segment in the template being reversed
    @property bool mate_is_reverse_strand()   const nothrow { return cast(bool)(flag & 0x20); }
    /// The first segment in the template
    @property bool is_first_of_pair()         const nothrow { return cast(bool)(flag & 0x40); }
    /// The last segment in the template
    @property bool is_second_of_pair()        const nothrow { return cast(bool)(flag & 0x80); }
    /// Secondary alignment
    @property bool is_secondary_alignment()   const nothrow { return cast(bool)(flag & 0x100); }
    /// Not passing quality controls
    @property bool failed_quality_control()   const nothrow { return cast(bool)(flag & 0x200); }
    /// PCR or optical duplicate
    @property bool is_duplicate()             const nothrow { return cast(bool)(flag & 0x400); }

    // flag setters
    @property void is_paired(bool b)                { _setFlag( 0, b); }
    @property void proper_pair(bool b)              { _setFlag( 1, b); }
    @property void is_unmapped(bool b)              { _setFlag( 2, b); }
    @property void mate_is_unmapped(bool b)         { _setFlag( 3, b); } 
    @property void is_reverse_strand(bool b)        { _setFlag( 4, b); } 
    @property void mate_is_reverse_strand(bool b)   { _setFlag( 5, b); } 
    @property void is_first_of_pair(bool b)         { _setFlag( 6, b); } 
    @property void is_second_of_pair(bool b)        { _setFlag( 7, b); } 
    @property void is_secondary_alignment(bool b)   { _setFlag( 8, b); } 
    @property void failed_quality_control(bool b)   { _setFlag( 9, b); } 
    @property void is_duplicate(bool b)             { _setFlag(10, b); } 

    @property string read_name() const nothrow {
        // notice -1: the string is zero-terminated, so we should strip that '\0'
        return cast(string)(_chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1]);
    }

    @property void read_name(string name) {
        enforce(name.length >= 1 && name.length <= 255, "name length must be in 1-255 range");
        _dup();
        utils.array.replaceSlice(_chunk, 
                     _chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1],
                     cast(ubyte[])name);
        _l_read_name = cast(ubyte)(name.length + 1);
    }

    /// List of CIGAR operations
    @property const(CigarOperation)[] cigar() const nothrow {
        return cast(const(CigarOperation)[])(_chunk[_cigar_offset .. _cigar_offset + 
                                             _n_cigar_op * CigarOperation.sizeof]);
    }

    @property void cigar(const(CigarOperation)[] c) {
        _dup();
        utils.array.replaceSlice(_chunk,
             _chunk[_cigar_offset .. _cigar_offset + _n_cigar_op * CigarOperation.sizeof],
             cast(ubyte[])c);

        _n_cigar_op = cast(ushort)(c.length);

        _recalculate_bin();
    }

    /// The number of reference bases covered
    int bases_covered() {

        if (this.is_unmapped) {
            return 0; // actually, valid alignments should have empty cigar string
        }

        int n = 0;
        foreach (c; cigar) {
            switch (c.operation) {
                case 'M':
                case '=':
                case 'X':
                case 'D':
                case 'N':
                    n += c.length;
                    break;
                default:
                    break;
            }
        }
        return n;
    }

    /// Human-readable representation of CIGAR string
    string cigarString() {
        char[] str;

        // guess size of resulting string
        str.reserve(_n_cigar_op * 3);

        foreach (cigar_op; cigar) {
            str ~= to!string(cigar_op.length);
            str ~= cigar_op.operation;
        }
        return cast(string)str;
    }

    /// Sequence data 
    @property const(ubyte)[] raw_sequence_data() const nothrow {
        return _chunk[_seq_offset .. _seq_offset + (_l_seq + 1) / 2];
    }

    /// Range of characters
    @property auto sequence() const {
        auto even = map!"a>>4"(raw_sequence_data);
        auto odd  = map!"a&15"(raw_sequence_data);
        return map!`"=ACMGRSVTWYHKDBN"[a]`(take(roundRobin(even, odd), sequence_length));
    }

    /// Set sequence. Must be of the same length as current sequence.
    @property void sequence(string seq) {

        enforce(seq.length == _l_seq, "Sequence must have the same length as current");

        _dup();

        /// Translates sequence character to its internal representation.
        static ubyte toHalfByte(char c) {
            /// The table is taken from samtools/bam_import.c
            static ubyte[256] table = [
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
                15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
                15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
                15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
                15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
                15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
            ];
            return table[c];
        }

        auto _raw_length = (seq.length + 1) / 2;
        // set sequence
        ubyte[] _seq = _chunk[_seq_offset .. _seq_offset + _raw_length];
        for (size_t i = 0; i < _raw_length; ++i) {
            _seq[i] = cast(ubyte)(toHalfByte(seq[2 * i]) << 4);

            if (seq.length > 2 * i + 1)
                _seq[i] |= cast(ubyte)(toHalfByte(seq[2 * i + 1]));
        }
    }

    /// Quality data
    @property const(ubyte)[] phred_base_quality() const nothrow {
        return _chunk[_qual_offset .. _qual_offset + _l_seq * char.sizeof];
    }

    /// Set quality data - array length must be of the same length as the sequence.
    @property void phred_base_quality(const(ubyte)[] quality) {
        enforce(quality.length == _l_seq, "Quality data must be of the same length as sequence");
        _dup();
        _chunk[_qual_offset .. _qual_offset + _l_seq] = quality;
    }

    /**
      Constructs the struct from memory chunk
      */
    this(ubyte[] chunk) {

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

        this._is_slice = true;

        _chunk = chunk;
        if (std.system.endian != Endian.littleEndian) {
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
            // (after them read_name follows which is string)      
            //                                                     
            switchEndianness(_chunk.ptr, 8 * uint.sizeof);

            // Then we need to switch endianness of CIGAR data:
            switchEndianness(_chunk.ptr + _cigar_offset, 
                             _n_cigar_op * uint.sizeof);

            // Dealing with tags is the responsibility of TagStorage.
        }

        // create from chunk of little-endian memory
        if (std.system.endian != Endian.littleEndian) {
            fixTagStorageByteOrder();
        }
    } 
   
    /// Construct alignment from basic information about it.
    ///
    /// Other fields can be set afterwards.
    this(string read_name,                          // info for developers:
         string sequence,                           // these 3 fields are needed
         in CigarOperation[] cigar)                 // to calculate size of _chunk
    {
        enforce(read_name.length < 256, "Too long read name, length must be <= 255");

        this._is_slice = false;

        this._chunk = new ubyte[8 * int.sizeof
                              + (read_name.length + 1) // tailing '\0'
                              + uint.sizeof * cigar.length
                              + ubyte.sizeof * ((sequence.length + 1) / 2)
                              + ubyte.sizeof * sequence.length];
        
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

        this.sequence = sequence;
    }

    /// Compare two alignments, including tags 
    /// (they must be in the same order for equality).
    bool opEquals(const ref Alignment other) const pure nothrow {
        // in D, array comparison compares elements, not pointers.
        return this._chunk == other._chunk;
    }

    /// Size of alignment when output to stream in BAM format.
    /// Includes block_size as well (see SAM/BAM specification)
    @property auto size_in_bytes() const {
        return int.sizeof + _chunk.length;
    }
   
    /// Write alignment to EndianStream, together with block_size
    /// and auxiliary data.
    void write(EndianStream stream) {
        stream.write(cast(int)(_chunk.length));
		stream.writeExact(_chunk.ptr, _tags_offset);
        writeTags(stream);
    }

private:

    ubyte[] _chunk; /// holds all the data, 
                    /// the access is organized via properties
                    /// (see below)

    bool _is_slice; /// indicates whether _chunk is a slice or an allocated array.

    /// Official field names from SAM/BAM specification.
    /// For internal use only

    @property  int _refID()      const nothrow { 
        return *(cast( int*)(_chunk.ptr + int.sizeof * 0)); 
    }

    @property  int _pos()        const nothrow { 
        return *(cast( int*)(_chunk.ptr + int.sizeof * 1)); 
    }

    @property uint _bin_mq_nl()  const nothrow { 
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

    /// Setters, also only for internal use
    @property void _refID(int n)       { *(cast( int*)(_chunk.ptr + int.sizeof * 0)) = n; }
    @property void _pos(int n)         { *(cast( int*)(_chunk.ptr + int.sizeof * 1)) = n; }
    @property void _bin_mq_nl(uint n)  { *(cast(uint*)(_chunk.ptr + int.sizeof * 2)) = n; }
    @property void _flag_nc(uint n)    { *(cast(uint*)(_chunk.ptr + int.sizeof * 3)) = n; }
    @property void _l_seq(int n)       { *(cast( int*)(_chunk.ptr + int.sizeof * 4)) = n; }
    @property void _next_refID(int n)  { *(cast( int*)(_chunk.ptr + int.sizeof * 5)) = n; }
    @property void _next_pos(int n)    { *(cast( int*)(_chunk.ptr + int.sizeof * 6)) = n; }
    @property void _tlen(int n)        { *(cast( int*)(_chunk.ptr + int.sizeof * 7)) = n; }

    /// Additional useful properties, also from SAM/BAM specification
    ///
    ///             The layout of bin_mq_nl and flag_nc is as follows
    ///                     (upper bits -------> lower bits):
    /// 
    /// bin_mq_nl [ { bin (16b) }  { mapping quality (8b) } { read name length (8b) } ]
    ///
    /// flag_nc   [ { flag (16b) } { n_cigar_op (16b) } ]
    ///
    @property ushort _bin()         const nothrow { 
        return _bin_mq_nl >> 16; 
    }
    @property  ubyte _mapq()        const nothrow { 
        return (_bin_mq_nl >> 8) & 0xFF; 
    }
    @property  ubyte _l_read_name() const nothrow { 
        return _bin_mq_nl & 0xFF; 
    }
    @property ushort _flag()        const nothrow { 
        return _flag_nc >> 16; 
    }
    @property ushort _n_cigar_op()  const nothrow { 
        return _flag_nc & 0xFFFF; 
    }
  
    /// Setters for those properties
    @property void _bin(ushort n)         { _bin_mq_nl = (_bin_mq_nl &  0xFFFF) | (n << 16); } 
    @property void _mapq(ubyte n)         { _bin_mq_nl = (_bin_mq_nl & ~0xFF00) | (n << 8); }
    @property void _l_read_name(ubyte n)  { _bin_mq_nl = (_bin_mq_nl & ~0xFF  ) | n; }
    @property void _flag(ushort n)        { _flag_nc   = (_flag_nc   &  0xFFFF) | (n << 16); }
    @property void _n_cigar_op(ushort n)  { _flag_nc   = (_flag_nc   & ~0xFFFF) | n; }

    /// Offsets of various arrays in bytes.
    /// Currently, are computed each time, so if speed will be an issue,
    /// they can be made fields instead of properties.
    @property size_t _read_name_offset() const nothrow { 
        return 8 * int.sizeof; 
    }

    @property size_t _cigar_offset()     const nothrow { 
        return _read_name_offset + _l_read_name * char.sizeof; 
    }

    @property size_t _seq_offset()       const nothrow { 
        return _cigar_offset + _n_cigar_op * uint.sizeof; 
    }

    @property size_t _qual_offset()      const nothrow { 
        return _seq_offset + (_l_seq + 1) / 2 * ubyte.sizeof; 
    }
    /// Offset of auxiliary data
    @property size_t _tags_offset()      const nothrow { 
        return _qual_offset + _l_seq * char.sizeof; 
    }

    /// Sets n-th flag bit to boolean value b.
    void _setFlag(int n, bool b) {
        assert(n < 16);
        // http://graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
        ushort mask = cast(ushort)(1 << n);
        _flag = (_flag & ~mask) | ((-cast(int)b) & mask);
    }

    /// If _chunk is still a slice, not an array, duplicate it.
    /// Used when some part of alignment record is modified by user.
    ///
    /// Basically, it's sort of copy-on-write: a lot of read-only alignments
    /// may point to the same location, but every modified one allocates its
    /// own chunk of memory.
    void _dup() {
        if (_is_slice) {
            _is_slice = false;
            _chunk = _chunk.dup;
        }
    }

    /// Calculates bin number.
    void _recalculate_bin() {
        _bin = reg2bin(position, position + bases_covered());
    }
}


///////////////////////////////////////////////////////////////////////////////
///
/// Lazy tag storage. 
///
///   Provides hash-like access (currently read-only) and opportunity to iterate
///   storage like an associative array.
///
///////////////////////////////////////////////////////////////////////////////
mixin template TagStorage() {

    ////////////////////////////////////////////////////////////////////////////
    ///
    /// Provides access to chunk of memory which contains tags
    /// This way, every time _tags_offset gets updated
    /// (due to update of cigar string of read name and memory move),
    /// the change is reflected automatically in tag storage.
    ///
    ////////////////////////////////////////////////////////////////////////////
    private @property const(ubyte)[] _tags_chunk() const  {
        return _chunk[_tags_offset .. $];
    }

    ////////////////////////////////////////////////////////////////////////////
    ///
    ///  Hash-like access to tags.
    ///
    ////////////////////////////////////////////////////////////////////////////
    Value opIndex(string key) {
        enforce(key.length == 2, "Key length must be 2");
        if (_tags_chunk.length < 4)
            return Value(null);
        
       size_t offset = 0;
       while (offset + 1 < _tags_chunk.length) {
           if (_tags_chunk[offset .. offset + 2] == key) {
               offset += 2;
               return readValue(offset);
           } else {
               offset += 2;
               skipValue(offset);
           }
       }
       return Value(null);
    }

    /// ditto
    void opIndexAssign(Value value, string key) {
        enforce(key.length == 2, "Key length must be 2");

        _dup();

        size_t offset = 0;
        while (offset + 1 < _tags_chunk.length) {
            if (_tags_chunk[offset .. offset + 2] == key) {
                replaceValueAt(offset + 2, value);
                return;
            } else {
                offset += 2;
                skipValue(offset);
            }
        }

        // append new tag to the end
        auto oldlen = _chunk.length;
        _chunk.length = _chunk.length + sizeInBytes(value) + 2 * char.sizeof;
        _chunk[oldlen .. oldlen + 2] = cast(ubyte[])key;
        emplaceValue(_chunk.ptr + oldlen + 2, value);
    }

    /// replace existing tag
    private void replaceValueAt(size_t offset, Value value) {
        // offset points to the beginning of the value
        auto begin = offset;
        skipValue(offset); // now offset is updated and points to the end
        auto end = offset;
        
        prepareSlice(_chunk, _tags_chunk[begin .. end], sizeInBytes(value));

        emplaceValue(_chunk.ptr + _tags_offset + begin, value);
    }

    /////////////////////////////////////////////////////////////////////////////
    ///  Provides opportunity to iterate over tags.
    /////////////////////////////////////////////////////////////////////////////
    int opApply(int delegate(ref string k, ref Value v) dg) {
        size_t offset = 0;
        while (offset + 1 < _tags_chunk.length) {
            auto key = cast(string)_tags_chunk[offset .. offset + 2];
            offset += 2;
            auto val = readValue(offset);
            auto res = dg(key, val);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

    /// Writes auxiliary data to output stream
    private void writeTags(Stream stream) {
		if (std.system.endian == Endian.littleEndian) {
			stream.writeExact(_tags_chunk.ptr, _tags_chunk.length);
		} else {
			fixTagStorageByteOrder();                                // FIXME: should modify on-the-fly
			stream.writeExact(_tags_chunk.ptr, _tags_chunk.length);  // during writing to the stream
			fixTagStorageByteOrder();                                
		}
    }

    ////////////////////////////////////////////////////////////////////////////
    ///  Reads value which starts from (_tags_chunk.ptr + offset) address,
    ///  and updates offset to the end of value.
    ////////////////////////////////////////////////////////////////////////////
    private Value readValue(ref size_t offset) {

        string readValueArrayTypeHelper() {
            char[] cases;
            foreach (c2t; ArrayElementTagValueTypes) {
                cases ~= 
                "case '"~c2t.ch~"':".dup~
                "  auto begin = offset;"~
                "  auto end = offset + length * "~c2t.ValueType.stringof~".sizeof;"~
                "  offset = end;"~ 
                "  return Value(cast("~c2t.ValueType.stringof~"[])(_tags_chunk[begin .. end].dup));";
                // TODO: copy-on-write in Value type (currently, dup is used)
            }
            return to!string("switch (elem_type) {" ~ cases ~
                   "  default: throw new UnknownTagTypeException(to!string(elem_type));"~
                   "}");
        }

        string readValuePrimitiveTypeHelper() {
            char[] cases;
            foreach (c2t; PrimitiveTagValueTypes) {
                cases ~= "case '"~c2t.ch~"':"~
                         "  auto p = _tags_chunk.ptr + offset;"~ 
                         "  auto value = *(cast("~c2t.ValueType.stringof~"*)p);"~
                         "  offset += value.sizeof;"~
                         "  return Value(value);".dup;
            }
            return to!string("switch (type) {" ~ cases ~
                   "  default: throw new UnknownTagTypeException(to!string(type));"~
                   "}");
        }

        char type = cast(char)_tags_chunk[offset++];
        if (type == 'Z' || type == 'H') {
            auto begin = offset;
            while (_tags_chunk[offset++] != 0) {}
            // return string with stripped '\0'
            auto v = Value(cast(string)_tags_chunk[begin .. offset - 1]);
            if (type == 'H') {
                v.setHexadecimalFlag();
            }
            return v;
        } else if (type == 'B') {
            char elem_type = cast(char)_tags_chunk[offset++];
            uint length = *(cast(uint*)(_tags_chunk.ptr + offset));
            offset += uint.sizeof;
            mixin(readValueArrayTypeHelper());
        } else {
            mixin(readValuePrimitiveTypeHelper());
        }
    }

    /**
      Increases offset so that it points to the next value.
    */
    private void skipValue(ref size_t offset) {
        char type = cast(char)_tags_chunk[offset++];
        if (type == 'Z' || type == 'H') {
            while (_tags_chunk[offset++] != 0) {}
        } else if (type == 'B') {
            char elem_type = cast(char)_tags_chunk[offset++];
            auto length = *(cast(uint*)(_tags_chunk.ptr + offset));
            offset += uint.sizeof + charToSizeof(elem_type) * length;
        } else {
            offset += charToSizeof(type);
        }
    }

    /**
      Intended to be used in constructor for initial endianness fixing
      in case the library is used on big-endian system.

      NOT TESTED AT ALL!!!
    */
    private void fixTagStorageByteOrder() {
        /* TODO: TEST ON BIG-ENDIAN SYSTEM!!! */
        const(ubyte)* p = _tags_chunk.ptr;
        const(ubyte)* end = _tags_chunk.ptr + _chunk.length;
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

    auto read = Alignment("readname", 
                          "AGCTGACTACGTAATAGCCCTA", 
                          [CigarOperation(22, 'M')]);
    assert(read.sequence_length == 22);
    assert(read.cigar.length == 1);
    assert(read.cigarString() == "22M");
    assert(read.read_name == "readname");
    assert(equal(read.sequence(), "AGCTGACTACGTAATAGCCCTA"));

    read.read_name = "anothername";
    assert(read.read_name == "anothername");
    assert(read.cigarString() == "22M");

    read.phred_base_quality = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                    13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
    assert(reduce!"a+b"(0, read.phred_base_quality) == 253);

    read.cigar = [CigarOperation(20, 'M'), CigarOperation(2, 'X')];
    assert(read.cigarString() == "20M2X");

    read.sequence = "AGCTGGCTACGTAATAGCCCTA";
    assert(equal(read.sequence(), "AGCTGGCTACGTAATAGCCCTA"));

    read["RG"] = Value(15);
    assert(to!int(read["RG"]) == 15);

    read["X1"] = Value([1, 2, 3, 4, 5]);
    assert(to!(int[])(read["X1"]) == [1, 2, 3, 4, 5]);

    read["RG"] = Value(5.6);
    assert(approxEqual(to!float(read["RG"]), 5.6));

    read["X1"] = Value(42);
    assert(to!int(read["X1"]) == 42);
}
