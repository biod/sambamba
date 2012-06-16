%%{
    machine sam_alignment;
 
    action update_sign { current_sign = fc == '-' ? -1 : 1; }
    action init_integer { int_value = 0; }
    action consume_next_digit { int_value *= 10; int_value += fc - '0'; }
    action take_sign_into_account { int_value *= current_sign; current_sign = 1; }

    sign = [\-+];

    uint = ([0-9]{1,18}) > init_integer $ consume_next_digit ;
    int = (sign >update_sign)? uint % take_sign_into_account ;

    action mark_float_start { float_beg = p - line.ptr; }
    action update_float_value { 
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }

    float = (sign? digit* '.'? digit+ ([eE] sign? digit+)?) > mark_float_start % update_float_value ;

    action qname_start { read_name_beg = p - line.ptr; }
    action qname_end { read_name_end = p - line.ptr; }

    qname =  '*' | ([!-?A-~]{1,255} > qname_start % qname_end) ;

    action set_flag { flag = to!ushort(int_value); }
    action set_pos { pos = to!uint(int_value); }
    action set_mapping_quality { mapping_quality = to!ubyte(int_value); }

    flag = uint % set_flag;
    rname = '*' | [!-()+-<>-~] [!-~]* ;
    pos = uint % set_pos;
    mapq = uint % set_mapping_quality;

    action cigar_set_op_length { cigar_op_len = to!uint(int_value); }
    action cigar_set_op_chr { cigar_op_chr = fc; }
    action cigar_put_operation { 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }

    cigar = '*' | (uint % cigar_set_op_length
                  [MIDNSHPX=] > cigar_set_op_chr % cigar_put_operation)+ ;

    action set_mate_pos { mate_pos = to!uint(int_value); }
    action set_template_length { template_length = to!int(int_value); }

    rnext = '*' | '=' | [!-()+-<>-~][!-~]* ;
    pnext = uint % set_mate_pos;
    tlen = int % set_template_length;

    action sequence_start { sequence_beg = p - line.ptr; }
    action sequence_end { sequence = line[sequence_beg .. p - line.ptr]; }

    seq = '*' | ([A-Za-z=.]+ > sequence_start % sequence_end) ;
    qual = [!-~]+ ;

    action create_alignment_struct {
        read = Alignment(line[read_name_beg .. read_name_end],
                         sequence,
                         cigar.data);
    }

    action allocate_quality_array {
        if (sequence.length > 1024) {
            qual_ptr = (new ubyte[sequence.length]).ptr;
        } else {
            qual_ptr = cast(ubyte*)alloca(sequence.length);
            if (!qual_ptr) {
                qual_ptr = (new ubyte[sequence.length]).ptr;
            }
        }
        qual_index = 0;
    }

    action convert_next_character_to_prob {
        qual_ptr[qual_index++] = cast(ubyte)(fc - 33);
    }

    action set_quality_data {
        read.phred_base_quality = qual_ptr[0 .. sequence.length];
    }

    mandatoryfields = qname '\t'
                       flag '\t'
                      rname '\t'
                        pos '\t'
                       mapq '\t'
                      cigar '\t'
                      rnext '\t'
                      pnext '\t'
                       tlen '\t'
                       (seq '\t' % create_alignment_struct)
                      (qual > allocate_quality_array
                            $ convert_next_character_to_prob
                            % set_quality_data);

    action set_charvalue { current_tagvalue = Value(fc); }
    action set_integervalue { 
        if (int_value < 0) {
            if (int_value >= byte.min) {
                current_tagvalue = Value(to!byte(int_value));
            } else if (int_value >= short.min) {
                current_tagvalue = Value(to!short(int_value));
            } else if (int_value >= int.min) {
                current_tagvalue = Value(to!int(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                current_tagvalue = Value(to!ubyte(int_value));
            } else if (int_value <= ushort.max) {
                current_tagvalue = Value(to!ushort(int_value));
            } else if (int_value <= uint.max) {
                current_tagvalue = Value(to!uint(int_value));
            } else {
                throw new Exception("integer out of range");
            }
        }
    }

    action start_tagvalue { tagvalue_beg = p - line.ptr; }

    action set_floatvalue { 
        current_tagvalue = Value(float_value);
    }

    action set_stringvalue { 
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]); 
    }

    action set_hexstringvalue {
        current_tagvalue = Value(line[tagvalue_beg .. p - line.ptr]);
        current_tagvalue.setHexadecimalFlag();
    }

    charvalue =  [!-~] > set_charvalue ;
    integervalue = int % set_integervalue;
    floatvalue = float % set_floatvalue ;

    action start_arrayvalue {
        // it might be not the best idea to use outbuffer;
        // the better idea might be two-pass approach
        // when first pass is for counting commas, and
        // the second is for filling allocated array
        outbuffer = new OutBuffer();        
        arraytype = fc;
    }

    action put_integer_to_array {
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': outbuffer.write(to!byte(int_value)); break;
            case 'C': outbuffer.write(to!ubyte(int_value)); break;
            case 's': outbuffer.write(to!short(int_value)); break;
            case 'S': outbuffer.write(to!ushort(int_value)); break;
            case 'i': outbuffer.write(to!int(int_value)); break;
            case 'I': outbuffer.write(to!uint(int_value)); break;
            default: assert(0);
        }
    }

    action put_float_to_array { 
        outbuffer.write(float_value); 
    }

    action set_arrayvalue {
        switch (arraytype) {
            case 'c': current_tagvalue = Value(cast(byte[])(outbuffer.toBytes())); break;
            case 'C': current_tagvalue = Value(cast(ubyte[])(outbuffer.toBytes())); break;
            case 's': current_tagvalue = Value(cast(short[])(outbuffer.toBytes())); break;
            case 'S': current_tagvalue = Value(cast(ushort[])(outbuffer.toBytes())); break;
            case 'i': current_tagvalue = Value(cast(int[])(outbuffer.toBytes())); break;
            case 'I': current_tagvalue = Value(cast(uint[])(outbuffer.toBytes())); break;
            case 'f': current_tagvalue = Value(cast(float[])(outbuffer.toBytes())); break;
            default: assert(0);
        }
    }

    stringvalue = [ !-~]+ > start_tagvalue % set_stringvalue ;
    hexstringvalue = xdigit+ > start_tagvalue % set_hexstringvalue ;
    integerarrayvalue = [cCsSiI] > start_arrayvalue (',' int % put_integer_to_array)+ % set_arrayvalue;
    floatarrayvalue = [f] > start_arrayvalue (',' float % put_float_to_array)+ % set_arrayvalue;
    arrayvalue = integerarrayvalue | floatarrayvalue ;

    tagvalue = ("A:" charvalue) | 
               ("i:" integervalue) | 
               ("f:" floatvalue) | 
               ("Z:" stringvalue) | 
               ("H:" hexstringvalue) |
               ("B:" arrayvalue) ;

    action tag_key_start { tag_key_beg = p - line.ptr; }
    action tag_key_end   { current_tag = line[tag_key_beg .. p - line.ptr]; }
    action append_tag_value { read.appendTag(current_tag, current_tagvalue); }

    tag = (alpha alnum) > tag_key_start % tag_key_end ;
    optionalfield = tag ':' tagvalue % append_tag_value ;
    optionalfields = optionalfield ('\t' optionalfield)* ;

    alignment := mandatoryfields ('\t' optionalfields)? ;

    write data; 
}%%

import alignment;
import tagvalue;

import std.array;
import std.conv;
import std.typecons;
import std.outbuffer;
import std.c.stdlib;

Alignment parseAlignmentLine(string line) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    byte current_sign = 1;

    size_t read_name_beg;
    size_t read_name_end;

    size_t sequence_beg;
    string sequence;

    uint cigar_op_len;
    char cigar_op_chr;

    size_t cigar_op_len_start;
    
    auto cigar = appender!(CigarOperation[])();

    long int_value; 
    float float_value;
    size_t float_beg;
    OutBuffer outbuffer;
    char arraytype;

    Alignment read;

    ushort flag;
    uint pos;
    uint mate_pos;
    ubyte mapping_quality; 
    int template_length;
    ubyte* qual_ptr;
    size_t qual_index;

    string current_tag;
    Value current_tagvalue;

    size_t tag_key_beg, tagvalue_beg;
    // TODO: RNAME, RNEXT

    %%write init;
    %%write exec;

    read.flag = flag;
    read.mapping_quality = mapping_quality;
    read.position = pos;
    read.template_length = template_length;
    read.next_pos = mate_pos;

    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto alignment = parseAlignmentLine(line);
    assert(alignment.read_name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60033);
    assert(alignment.mapping_quality == 25);
    assert(alignment.next_pos == 60033);

    import std.stdio;
    import sam.serialize;

    assert(to!char(alignment["XT"]) == 'U');
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
}
