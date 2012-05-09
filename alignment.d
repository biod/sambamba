module alignment;

import std.stream;
import std.algorithm;
import std.conv;
import std.system;
import std.variant;

/**
  Represents tag value.

  Algebraic type, includes all types mentioned in SAM/BAM specification.

  A -> char
  c -> byte
  C -> ubyte
  s -> short
  S -> ushort
  i -> int
  I -> uint
  f -> float
  H -> string
  Z -> string
  B -> array of [cCsSiIf]
  */
alias Algebraic!(char, byte, ubyte, short, ushort, int, uint, float,
                 string, byte[], ubyte[], short[], ushort[],
                 int[], uint[], float[]) Value;

/**
  Thrown in case of unrecognized tag type
 */
class UnknownTagTypeException : Exception {
    this(string msg) { super(msg); }
}

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
        if ((raw & 0xF) >= 8) {
            return '\0';
        }

        return "MIDNSHP=X"[raw & 0xF];
    }
}

private {
    struct when(char c, T) {
        enum ch = c;
        alias T ValueType;
    }

    Value readPrimitive(W...)(ref Stream stream, char type) {
        string readPrimitiveHelper() {
            char[] cases;
            foreach (w; W) {
                static if (is(w.ValueType == string)) {
                    // read zero-terminated string
                    /* TODO: optimize, avoid too much memory allocations */
                    cases ~= "case '".dup~w.ch~"':"~
                             "    char[] s;"~
                             "    char c;"~
                             "    do { stream.read(c); s ~= c; } while (c != 0);"~
                             "    return Value(cast(string)s[0..$-1]);";
                } else {
                    cases ~= "case '".dup~w.ch~"':"~
                             "    "~w.ValueType.stringof~" val;"~
                             "    stream.read(val);"~
                             "    return Value(val);";
                }
            }
            return to!string("switch (type) {"~ 
                        cases~
                   "    default: throw new UnknownTagTypeException(to!string(type));"~
                   "}");
        }
        mixin(readPrimitiveHelper());
    }

    Value readArray(W...)(ref Stream stream, char elem_type) {
        string readArrayHelper() {
            char[] cases;
            foreach (w; W) {
                auto w_str = w.ValueType.stringof;
                cases ~= "case '"~w.ch~"':"~
                         "    auto val = new "~w_str~"[length];"~
                         "    foreach (i; 0 .. length) {"~
                         "        stream.read(val[i]);"~
                         "    }"~
                         "    return Value(val);";
            }
            return to!string("switch (elem_type) {"~
                        cases~
                   "    default: throw new UnknownTagTypeException(elem_type ~ \"[]\");"~
                   "}");
        }
        uint length;
        stream.read(length);
        mixin(readArrayHelper());
    }
}

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

    CigarOperation[] cigar = void;

    ubyte[] seq = void;
    char[] qual = void;

    Value[string] tags;

    static Alignment parse(char[] chunk) {
        scope Stream memory_stream = new MemoryStream(chunk);
        scope Stream stream = new EndianStream(memory_stream,
                                         Endian.littleEndian);

        Alignment a;

        stream.read(a.ref_id);
        stream.read(a.position);

        uint bin_mq_nl = void;
        stream.read(bin_mq_nl);
        a.bin = bin_mq_nl >> 16;
        a.mapping_quality = (bin_mq_nl >> 8) & 0xFF;
        
        ubyte l_read_name = bin_mq_nl & 0xFF;

        uint flag_nc = void;
        stream.read(flag_nc);
        a.flag = flag_nc >> 16;

        ushort n_cigar_op = flag_nc & 0xFFFF;
        
        stream.read(a.sequence_length);
        stream.read(a.next_ref_id);
        stream.read(a.next_pos);
        stream.read(a.template_length);

        a.read_name = stream.readString(l_read_name)[0..$-1].idup;

        a.cigar = new CigarOperation[n_cigar_op];
        foreach (i; 0 .. n_cigar_op) {
            stream.read(a.cigar[i].raw);
        }
        
        a.seq = cast(ubyte[])stream.readString((a.sequence_length + 1) / 2);
        a.qual = stream.readString(a.sequence_length);

        // read auxiliary data
        /* TODO:
           store the whole chunk of memory in the structure
           pros: 1) string slicing instead of copying => less memory allocations
                 2) on little endian architectures, arrays can also be sliced
           */
        while (!stream.eof()) {
            char[] tag = stream.readString(2);
            char type = void;
            stream.read(type);

            if (type != 'B') {
                Value val = readPrimitive!(
                    when!('A', char),
                    when!('c', byte),  when!('C', ubyte),
                    when!('s', short), when!('S', ushort),
                    when!('i', int),   when!('I', uint),
                    when!('f', float),
                    when!('Z', string),
                    when!('H', string))(stream, type);
                a.tags[tag.idup] = val;
            } else {
                // array
                char elem_type = void;
                stream.read(elem_type);
                Value val = readArray!(
                    when!('c', byte),  when!('C', ubyte),
                    when!('s', short), when!('S', ushort),
                    when!('i', int),   when!('I', uint),
                    when!('f', float))(stream, elem_type);
                a.tags[tag.idup] = val;
            }
        } // end of alignment block
        return a;
    } // end of parse
}

auto alignmentRange(ref Stream stream) {
    return map!(Alignment.parse)(unparsedAlignments(stream));
}
