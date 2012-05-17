module alignment;
import tagvalue;
import tagstorage;

import std.stream;
import std.algorithm;
import std.system;

import utils.switchendianness;

/**
  Range for iterating over unparsed alignments
 */
class UnparsedAlignmentRange {

    this(ref Stream stream) {
        _stream = stream;
        readNext();
    }

    bool empty() @property {
        return _empty;
    }

    /**
        Returns: alignment block except the first 4 bytes (block_size)
     */
    ubyte[] front() @property {
        return _current_record;
    }

    void popFront() {
        readNext();
    }

private:
    Stream _stream;
    ubyte[] _current_record;
    bool _empty = false;

    /**
      Reads next alignment block from stream.
     */
    void readNext() {
        if (_stream.eof()) {
            _empty = true;
            return;
        }

        int block_size = void;
        _stream.read(block_size);
        _current_record = new ubyte[block_size];
        _stream.readExact(_current_record.ptr, block_size);
    }

}

/**
    Returns: range for iterating over alignment blocks
 */
auto unparsedAlignments(ref Stream stream) {
    return new UnparsedAlignmentRange(stream);
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

    string read_name() @property {
        // notice -1: the string is zero-terminated, so we should strip that '\0'
        return cast(string)(_chunk[_read_name_offset .. _read_name_offset + _l_read_name - 1]);
    }

    CigarOperation[] cigar() @property {
        return cast(CigarOperation[])(_chunk[_cigar_offset .. _cigar_offset + 
                                             _n_cigar_op * CigarOperation.sizeof]);
    }

    /// Human-readable representation of CIGAR string
    ///
    /// Evaluated lazily from raw data, therefore NOT
    /// marked as @property.
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
    /// Evaluated lazily, so not @property
    string sequence() {
        immutable string chars = "=ACMGRSVTWYHKDBN";
        char[] s = new char[sequence_length];
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
    char[] phred_base_quality() @property {
        return cast(char[])_chunk[_qual_offset .. _qual_offset + _l_seq * char.sizeof];
    }

    TagStorage tags = void;

    /**
      Constructs the struct from memory chunk
      */
    this(ubyte[] chunk) {

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

            // Dealing with tags is not our responsibility.
            // It goes to TagStorage interface implementations
            // which are given a reference to the chunk of memory
            // and are allowed to do whatever they find appropriate.
        }

        /// currently, LazyTagStorage is the best choice
        this.tags = new LazyTagStorage(cast(ubyte[])_chunk[_tags_offset .. $]);
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
  ulong _read_name_offset() @property { return 8 * int.sizeof; }
  ulong _cigar_offset()     @property { return _read_name_offset + _l_read_name * char.sizeof; }
  ulong _seq_offset()       @property { return _cigar_offset + _n_cigar_op * uint.sizeof; }
  ulong _qual_offset()      @property { return _seq_offset + (_l_seq + 1) / 2 * ubyte.sizeof; }

    /// Offset of auxiliary data
  ulong _tags_offset()      @property { return _qual_offset + _l_seq * char.sizeof; }
}

Alignment parseAlignment(ubyte[] chunk) {
    return Alignment(chunk);
}

import std.parallelism;

auto alignmentRange(ref Stream stream, TaskPool task_pool) {
    version(serial) {
        return map!parseAlignment(unparsedAlignments(stream));
    } else {
        /* TODO: tweak granularity */
        return map!parseAlignment(unparsedAlignments(stream));
//        return task_pool.map!parseAlignment(unparsedAlignments(stream), 500);
    }
}
