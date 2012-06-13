%%{
    machine sam_alignment;
 
    action update_sign { current_sign = fc == '-' ? -1 : 1; }
    action init_integer { int_value = 0; }
    action consume_next_digit { int_value *= 10; int_value += fc - '0'; }
    action take_sign_into_account { int_value *= current_sign; current_sign = 1; }

    sign = [\-+];

    uint = ([0-9]{1,9}) > init_integer $ consume_next_digit ;
    int = (sign >update_sign)? uint % take_sign_into_account ;

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

    charvalue =  [!-~] ;
    integervalue = int ;
    floatvalue = sign? digit* '.'? digit+ ([eE] sign? digit+)? ; 
    stringvalue = [ !-~]+ ;
    hexstringvalue = xdigit+ ;
    integerarrayvalue = [cCsSiI] (',' integervalue)+ ;
    floatarrayvalue = [f] (',' floatvalue)+ ;
    arrayvalue = integerarrayvalue | floatarrayvalue ;

    tagvalue = ("A:" charvalue) | 
               ("i:" integervalue) | 
               ("f:" floatvalue) | 
               ("Z:" stringvalue) | 
               ("H:" hexstringvalue) |
               ("B:" arrayvalue) ;

    tag = alpha alnum ; 
    optionalfield = tag ':' tagvalue ;
    optionalfields = optionalfield ('\t' optionalfield)* ;

    alignment := mandatoryfields ('\t' optionalfields)? ;

    write data; 
}%%

import alignment;
import std.array;
import std.conv;
import std.c.stdlib;

Alignment parseAlignmentLine(string line) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = null;
    int cs;

    byte current_sign = 1;

    size_t read_name_beg = 0;
    size_t read_name_end = 0;

    size_t sequence_beg = 0;
    string sequence;

    uint cigar_op_len;
    char cigar_op_chr;

    size_t cigar_op_len_start;
    
    auto cigar = appender!(CigarOperation[])();

    int int_value; 
    Alignment read;

    ushort flag;
    uint pos;
    uint mate_pos;
    ubyte mapping_quality; 
    int template_length;
    ubyte* qual_ptr;
    size_t qual_index;
   
    // TODO: RNAME, RNEXT, tags

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

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35 MD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";

	auto alignment = parseAlignmentLine(line);
    assert(alignment.read_name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60033);
    assert(alignment.mapping_quality == 25);
    assert(alignment.next_pos == 60033);
}
