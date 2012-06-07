module alignment;

import tagvalue;
import tagstorage;

private import bai.bin;

import std.algorithm;
import std.range;
import std.conv;
import std.exception;
import std.system;

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
    uint length() @property const {
        return raw >> 4;
    }
  
    /// CIGAR operation as one of MIDNSHP=X
    char operation() @property const {
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

    // /////////////////////////////////////////////////////////////////////////
    //
    // Layout of Alignment in memory:
    //
    // TagStorage      <- tags
    //     ubyte[]     <- tags/_chunk: slice of _chunk (see below) or an array
    //         size_t                  containing tags in binary form, as in BAM.
    //         ubyte*
    // ubyte[]         <- _chunk: binary representation of alignment, as in BAM.
    //     size_t
    //     ubyte*
    // bool            <- _is_slice: indicates whether _chunk is a slice or not
    //
    // //////////////////////////////////////////////////////////////////////////

    @property    int ref_id()           const { return _refID; }
    @property    int position()         const { return _pos; }
    @property ushort bin()              const { return _bin; }
    @property  ubyte mapping_quality()  const { return _mapq; }
    @property ushort flag()             const { return _flag; }
    @property    int sequence_length()  const { return _l_seq; }
    @property    int next_ref_id()      const { return _next_refID; }
    @property    int next_pos()         const { return _next_pos; }
    @property    int template_length()  const { return _tlen; }

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
    @property bool is_paired()                const { return cast(bool)(flag & 0x1); }
    /// Each segment properly aligned according to the aligner
    @property bool proper_pair()              const { return cast(bool)(flag & 0x2); }
    /// Segment unmapped
    @property bool is_unmapped()              const { return cast(bool)(flag & 0x4); }
    /// Next segment in the template unmapped
    @property bool mate_is_unmapped()         const { return cast(bool)(flag & 0x8); }
    /// Sequence being reverse complemented
    @property bool is_reverse_strand()        const { return cast(bool)(flag & 0x10); }
    /// Sequence of the next segment in the template being reversed
    @property bool mate_is_reverse_strand()   const { return cast(bool)(flag & 0x20); }
    /// The first segment in the template
    @property bool is_first_of_pair()         const { return cast(bool)(flag & 0x40); }
    /// The last segment in the template
    @property bool is_second_of_pair()        const { return cast(bool)(flag & 0x80); }
    /// Secondary alignment
    @property bool is_secondary_alignment()   const { return cast(bool)(flag & 0x100); }
    /// Not passing quality controls
    @property bool failed_quality_control()   const { return cast(bool)(flag & 0x200); }
    /// PCR or optical duplicate
    @property bool is_duplicate()             const { return cast(bool)(flag & 0x400); }

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

    @property string read_name() const {
        // notice -1: the string is zero-terminated, so we should strip that '\0'
        return cast(string)(_chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1]);
    }

    @property const(CigarOperation)[] cigar() const {
        return cast(const(CigarOperation)[])(_chunk[_cigar_offset .. _cigar_offset + 
                                             _n_cigar_op * CigarOperation.sizeof]);
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
    @property const(ubyte)[] raw_sequence_data() const {
        return _chunk[_seq_offset .. _seq_offset + (_l_seq + 1) / 2];
    }

    /// Range of characters
    auto sequence() const {
        auto even = map!"a>>4"(raw_sequence_data);
        auto odd  = map!"a&15"(raw_sequence_data);
        return map!`"=ACMGRSVTWYHKDBN"[a]`(take(roundRobin(even, odd), sequence_length));
    }

    /// Quality data
    @property const(ubyte)[] phred_base_quality() const {
        return _chunk[_qual_offset .. _qual_offset + _l_seq * char.sizeof];
    }

    TagStorage tags = void;

    /**
      Constructs the struct from memory chunk
      */
    this(ubyte[] chunk) {

        // TODO: switch endianness lazily as well?

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

        this.tags = TagStorage(_chunk[_tags_offset .. $]);
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
        
        this.tags = TagStorage(_chunk[$ .. $]);

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

        // set sequence
        auto _raw_length = (sequence.length + 1) / 2;
        ubyte[] _seq = _chunk[_seq_offset .. _seq_offset + _raw_length];
        for (size_t i = 0; i < _raw_length; ++i) {
            _seq[i] = cast(ubyte)(toHalfByte(sequence[2 * i]) << 4);

            if (sequence.length > 2 * i + 1)
                _seq[i] |= cast(ubyte)(toHalfByte(sequence[2 * i + 1]));
        }
    }

    bool opEquals(const ref Alignment other) const pure nothrow {
        /// Notice that in D, array comparison compares elements, not pointers.
        return this._chunk == other._chunk && this.tags == other.tags;
    }

private:

    ubyte[] _chunk; /// holds all the data, 
                    /// the access is organized via properties
                    /// (see below)

    bool _is_slice; /// indicates whether _chunk is a slice or an allocated array.

    /// Official field names from SAM/BAM specification.
    /// For internal use only

    @property  int _refID()      const { return *(cast( int*)(_chunk.ptr + int.sizeof * 0)); }
    @property  int _pos()        const { return *(cast( int*)(_chunk.ptr + int.sizeof * 1)); }
    @property uint _bin_mq_nl()  const { return *(cast(uint*)(_chunk.ptr + int.sizeof * 2)); }
    @property uint _flag_nc()    const { return *(cast(uint*)(_chunk.ptr + int.sizeof * 3)); }
    @property  int _l_seq()      const { return *(cast( int*)(_chunk.ptr + int.sizeof * 4)); }
    @property  int _next_refID() const { return *(cast( int*)(_chunk.ptr + int.sizeof * 5)); }
    @property  int _next_pos()   const { return *(cast( int*)(_chunk.ptr + int.sizeof * 6)); }
    @property  int _tlen()       const { return *(cast( int*)(_chunk.ptr + int.sizeof * 7)); }

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
    @property ushort _bin()         const { return _bin_mq_nl >> 16; }
    @property  ubyte _mapq()        const { return (_bin_mq_nl >> 8) & 0xFF; }
    @property  ubyte _l_read_name() const { return _bin_mq_nl & 0xFF; }
    @property ushort _flag()        const { return _flag_nc >> 16; }
    @property ushort _n_cigar_op()  const { return _flag_nc & 0xFFFF; }
  
    /// Setters for those properties
    @property void _bin(ushort n)         { _bin_mq_nl = (_bin_mq_nl &  0xFFFF) | (n << 16); } 
    @property void _mapq(ubyte n)         { _bin_mq_nl = (_bin_mq_nl & ~0xFF00) | (n << 8); }
    @property void _l_read_name(ubyte n)  { _bin_mq_nl = (_bin_mq_nl & ~0xFF  ) | n; }
    @property void _flag(ushort n)        { _flag_nc   = (_flag_nc   &  0xFFFF) | (n << 16); }
    @property void _n_cigar_op(ushort n)  { _flag_nc   = (_flag_nc   & ~0xFFFF) | n; }

    /// Offsets of various arrays in bytes.
    /// Currently, are computed each time, so if speed will be an issue,
    /// they can be made fields instead of properties.
    @property size_t _read_name_offset() const { return 8 * int.sizeof; }
    @property size_t _cigar_offset()     const { return _read_name_offset + _l_read_name * char.sizeof; }
    @property size_t _seq_offset()       const { return _cigar_offset + _n_cigar_op * uint.sizeof; }
    @property size_t _qual_offset()      const { return _seq_offset + (_l_seq + 1) / 2 * ubyte.sizeof; }

    /// Offset of auxiliary data
    @property size_t _tags_offset()      const { return _qual_offset + _l_seq * char.sizeof; }

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
    /// may copy to the same location, but every modified one allocates its
    /// own chunk of memory.
    ///
    /// Atomic compare-and-swap is used so multiple threads can modify fields.
    void _dup() {
        if (_is_slice) {
            _is_slice = false;
            _chunk = _chunk.dup;
            tags = TagStorage(_chunk[_tags_offset .. $]);
        }
    }

    /// Calculates bin number.
    void _recalculate_bin() {
        _bin = reg2bin(position, position + bases_covered());
    }
}

unittest {
    import std.algorithm;
    import std.stdio;
    auto read = Alignment("readname", 
                          "AGCTGACTACGTAATAGCCCTA", 
                          [CigarOperation(22, 'M')]);
    assert(read.sequence_length == 22);
    assert(read.cigar.length == 1);
    assert(read.cigarString() == "22M");
    assert(read.read_name == "readname");
    assert(equal(read.sequence(), "AGCTGACTACGTAATAGCCCTA"));
}
