module alignment;

import tagvalue;
import tagstorage;
import chunkinputstream;
import virtualoffset;

import std.stream;
import std.algorithm;
import std.system;
import std.array : uninitializedArray;
import std.conv;

import utils.switchendianness;

struct RawAlignmentBlock {
    VirtualOffset virtual_offset;
    ubyte[] raw_alignment_data;
}

/// Whether to provide offsets or not
/// when iterating raw alignment blocks
private enum IteratePolicy {
    withOffsets,
    withoutOffsets
}

/**
  Range for iterating over unparsed alignments,
  together with their virtual offsets in the BAM file.
 */
private class UnparsedAlignmentRange(IteratePolicy policy) 
{ 

    /// Create new range from IChunkInputStream.
    this(ref IChunkInputStream stream) {
        _stream = stream;
        _endian_stream = new EndianStream(_stream, Endian.littleEndian);
        readNext();
    }

    bool empty() @property {
        return _empty;
    }

static if (policy == IteratePolicy.withOffsets) {
    /**
        Returns: tuple of 1) virtual offset of alignment block beginning,
                          2) alignment block except the first 4 bytes (block_size)

        See SAM/BAM specification for details about the structure of alignment blocks.
     */
    RawAlignmentBlock front() @property {
        return RawAlignmentBlock(_current_voffset, _current_record);
    }
} else {
    /**
        Returns: alignment block except the first 4 bytes (block_size)
     */
    ubyte[] front() @property {
        return _current_record;
    }
}

    void popFront() {
        readNext();
    }

private:
    IChunkInputStream _stream;
    EndianStream _endian_stream;

    ubyte[] _current_record;
    bool _empty = false;

static if (policy == IteratePolicy.withOffsets) {
    VirtualOffset _current_voffset;
}

    /**
      Reads next alignment block from stream.
     */
    void readNext() {
        if (_stream.eof()) {
            _empty = true;
            return;
        }
        
static if (policy == IteratePolicy.withOffsets) {
        _current_voffset = _stream.virtualTell();
}

        int block_size = void;
        _endian_stream.read(block_size);
        _current_record = uninitializedArray!(ubyte[])(block_size);
        _endian_stream.readExact(_current_record.ptr, block_size);
    }

}

/**
    Returns: range for iterating over alignment blocks
 */
auto unparsedAlignments(alias Policy)(ref IChunkInputStream stream) {
    return new UnparsedAlignmentRange!(Policy)(stream);
}

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

    uint raw; /// raw data from BAM

    /// operation length
    uint length() @property {
        return raw >> 4;
    }
  
    /* FIXME: better name? */
    /// CIGAR operation as one of MIDNSHP=X
    char operation() @property {
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

    int ref_id() @property {
        return _refID;
    }

    int position() @property {
        return _pos;
    }

    ushort bin() @property {
        return _bin;
    }

    ubyte mapping_quality() @property {
        return (_bin_mq_nl >> 8) & 0xFF;
    }

    ushort flag() @property { 
        return _flag_nc >> 16;
    }

    int sequence_length() @property {
        return _l_seq;
    }

    int next_ref_id() @property {
        return _next_refID;
    }

    int next_pos() @property {
        return _next_pos;
    }

    int template_length() @property {
        return _tlen;
    }
  
    /// Template having multiple segments in sequencing
    bool is_paired()                @property { return cast(bool)(flag & 0x1); }
    /// Each segment properly aligned according to the aligner
    bool proper_pair()              @property { return cast(bool)(flag & 0x2); }
    /// Segment unmapped
    bool is_unmapped()              @property { return cast(bool)(flag & 0x4); }
    /// Next segment in the template unmapped
    bool mate_is_unmapped()         @property { return cast(bool)(flag & 0x8); }
    /// Sequence being reverse complemented
    bool is_reverse_strand()        @property { return cast(bool)(flag & 0x10); }
    /// Sequence of the next segment in the template being reversed
    bool mate_is_reverse_strand()   @property { return cast(bool)(flag & 0x20); }
    /// The first segment in the template
    bool is_first_of_pair()         @property { return cast(bool)(flag & 0x40); }
    /// The last segment in the template
    bool is_second_of_pair()        @property { return cast(bool)(flag & 0x80); }
    /// Secondary alignment
    bool is_secondary_alignment()   @property { return cast(bool)(flag & 0x100); }
    /// Not passing quality controls
    bool failed_quality_control()   @property { return cast(bool)(flag & 0x200); }
    /// PCR or optical duplicate
    bool is_duplicate()             @property { return cast(bool)(flag & 0x400); }

    string read_name() @property {
        // notice -1: the string is zero-terminated, so we should strip that '\0'
        return cast(string)(_chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1]);
    }

    CigarOperation[] cigar() @property {
        return cast(CigarOperation[])(_chunk[_cigar_offset .. _cigar_offset + 
                                             _n_cigar_op * CigarOperation.sizeof]);
    }

    /// The number of reference bases covered
    int bases_covered() {
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
    string cigar_string() {
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
    ubyte[] raw_sequence_data() @property {
        return _chunk[_seq_offset .. _seq_offset + (_l_seq + 1) / 2];
    }

    /// String representation of raw_sequence_data.
    string sequence() {
        immutable string chars = "=ACMGRSVTWYHKDBN";
        char[] s = uninitializedArray!(char[])(sequence_length);
        for (auto i = 0; i < sequence_length; i++) {
            auto j = i / 2;
            auto b = raw_sequence_data[j];
            if (i % 2 == 0) {
                s[i] = chars[b >> 4]; 
            } else {
                s[i] = chars[b & 0xF];
            }
        }
        return cast(string)s;
    }

    /// Quality data
    ubyte[] phred_base_quality() @property {
        return _chunk[_qual_offset .. _qual_offset + _l_seq * char.sizeof];
    }

    TagStorage tags = void;

    /**
      Constructs the struct from memory chunk
      */
    this(ubyte[] chunk) {

        // TODO: switch endianness lazily as well?

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

        this.tags = TagStorage(cast(ubyte[])_chunk[_tags_offset .. $]);
    } 

    bool opEquals(const ref Alignment other) const pure nothrow {
        return this._chunk == other._chunk && this.tags == other.tags;
    }

private:

    ubyte[] _chunk; /// holds all the data, 
                    /// the access is organized via properties
                    /// (see below)

    /// Official field names from SAM/BAM specification.
    /// For internal use only

    int _refID()            @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 0)); }
    int _pos()              @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 1)); }
   uint _bin_mq_nl()        @property { return *(cast(uint*)(_chunk.ptr + int.sizeof * 2)); }
   uint _flag_nc()          @property { return *(cast(uint*)(_chunk.ptr + int.sizeof * 3)); }
    int _l_seq()            @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 4)); }
    int _next_refID()       @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 5)); }
    int _next_pos()         @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 6)); }
    int _tlen()             @property { return *(cast( int*)(_chunk.ptr + int.sizeof * 7)); }

    /// Additional useful properties, also from SAM/BAM specification
 ushort _bin()              @property { return _bin_mq_nl >> 16; }
  ubyte _l_read_name()      @property { return _bin_mq_nl & 0xFF; }
 ushort _n_cigar_op()       @property { return _flag_nc & 0xFFFF; }
   
    /// Offsets of various arrays in bytes.
    /// Currently, are computed each time, so if speed will be an issue,
    /// they can be made fields instead of properties.
 size_t _read_name_offset() @property { return 8 * int.sizeof; }
 size_t _cigar_offset()     @property { return _read_name_offset + _l_read_name * char.sizeof; }
 size_t _seq_offset()       @property { return _cigar_offset + _n_cigar_op * uint.sizeof; }
 size_t _qual_offset()      @property { return _seq_offset + (_l_seq + 1) / 2 * ubyte.sizeof; }

    /// Offset of auxiliary data
 size_t _tags_offset()      @property { return _qual_offset + _l_seq * char.sizeof; }
}

/// Returns: an alignment constructed out of given chunk of memory.
///
/// No parsing occurs, it is done lazily.
Alignment makeAlignment(ubyte[] chunk) {
    return Alignment(chunk);
}

/// Tuple of virtual offset of the alignment, and the alignment itself.
struct AlignmentBlock {
    VirtualOffset virtual_offset;
    Alignment alignment;
}

/// Returns: lazy range of Alignment structs constructed from a given stream.
auto alignmentRange(ref IChunkInputStream stream) {

    return map!makeAlignment(unparsedAlignments!(IteratePolicy.withoutOffsets)(stream));
}

/// Returns: lazy range of AlignmentBlock structs constructed from a given stream.
///
/// Provides virtual offsets together with alignments and thus
/// Useful for random access operations.
auto alignmentRangeWithOffsets(ref IChunkInputStream stream) {
    static AlignmentBlock makeAlignmentBlock(RawAlignmentBlock block) {
        return AlignmentBlock(block.virtual_offset,
                              makeAlignment(block.raw_alignment_data));
    }

    return map!makeAlignmentBlock(unparsedAlignments!(IteratePolicy.withOffsets)(stream));
}
