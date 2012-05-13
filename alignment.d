module alignment;
import tagvalue;
import tagstorage;

import std.stream;
import std.algorithm;
import std.system;

private {

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
        char[] front() @property {
            return _current_record;
        }

        void popFront() {
            readNext();
        }

    private:
        Stream _stream;
        char[] _current_record;
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
            _current_record = _stream.readString(block_size);
        }
    }

    /**
        Returns: range for iterating over alignment blocks
     */
    auto unparsedAlignments(ref Stream stream) {
        return new UnparsedAlignmentRange(stream);
    }
}

/**
  Represents single CIGAR operation
 */
struct CigarOperation {
    uint raw; /// raw data from BAM

    /// operation length
    uint length() @property {
        return raw >> 4;
    }
   
    /// CIGAR operation as one of MIDNSHP=X
    char operation() @property {
        if ((raw & 0xF) > 8) {
            return '\0';
        }

        return "MIDNSHP=X"[raw & 0xF];
    }
}

/** 
  Represents single read
*/
struct Alignment {
    
    int ref_id = void;
    int position = void;
    ushort bin = void;
    ubyte mapping_quality = void;
    ushort flag = void;
    int sequence_length = void;

    int next_ref_id = void;
    int next_pos = void;
    int template_length = void;

    string read_name = void;

    CigarOperation[] cigar = void; /// CIGAR operations

    /// human-readable representation of CIGAR string
    string cigar_string() @property {
        char[] str;
        foreach (cigar_op; cigar) {
            str ~= to!string(cigar_op.length);
            str ~= cigar_op.operation;
        }
        return cast(string)str;
    }

    ubyte[] seq = void; /// sequence data

    /// string representation of seq
    string sequence() @property {
        immutable string chars = "=ACMGRSVTWYHKDBN";
        char[] s = new char[sequence_length];
        for (auto i = 0; i < sequence_length; i++) {
            auto j = i / 2;
            auto b = seq[j];
            if (i % 2 == 0) {
                s[i] = chars[b >> 4]; 
            } else {
                s[i] = chars[b & 0xF];
            }
        }
        return cast(string)s;
    }

    char[] qual = void; /// quality data

    TagStorage tags = void;

    /**
      Constructs the struct from memory chunk
      */
    this(char[] chunk) {
        scope Stream memory_stream = new MemoryStream(chunk);
        scope Stream stream = new EndianStream(memory_stream,
                                         Endian.littleEndian);

        stream.read(this.ref_id);
        stream.read(this.position);

        uint bin_mq_nl = void;
        stream.read(bin_mq_nl);
        this.bin = bin_mq_nl >> 16;
        this.mapping_quality = (bin_mq_nl >> 8) & 0xFF;
        
        ubyte l_read_name = bin_mq_nl & 0xFF;

        uint flag_nc = void;
        stream.read(flag_nc);
        this.flag = flag_nc >> 16;

        ushort n_cigar_op = flag_nc & 0xFFFF;
        
        stream.read(this.sequence_length);
        stream.read(this.next_ref_id);
        stream.read(this.next_pos);
        stream.read(this.template_length);

        this.read_name = stream.readString(l_read_name)[0..$-1].idup;

        this.cigar = new CigarOperation[n_cigar_op];
        foreach (i; 0 .. n_cigar_op) {
            stream.read(this.cigar[i].raw);
        }
       
        this.seq = cast(ubyte[])stream.readString((this.sequence_length + 1) / 2);
        this.qual = stream.readString(this.sequence_length);

        /* FIXME: don't hardcode storage type */
        this.tags = new LazyTagStorage(cast(ubyte[])chunk[cast(uint)stream.position .. $]);
    } 
}

Alignment parseAlignment(char[] chunk) {
    return Alignment(chunk);
}

import std.parallelism;

auto alignmentRange(ref Stream stream, TaskPool task_pool) {
    version(serial) {
        return map!parseAlignment(unparsedAlignments(stream));
    } else {
        /* TODO: tweak granularity */
        return task_pool.map!parseAlignment(unparsedAlignments(stream), 500);
    }
}
