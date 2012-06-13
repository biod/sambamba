%%{
    machine sam_alignment;
  
    sign = [\-+] ;

    action qname_start { read_name_beg = p - line.ptr; }
    action qname_end { read_name_end = p - line.ptr; }

    qname =  '*' | ([!-?A-~]{1,255} > qname_start % qname_end) ;
    flag = [0-9]{1,5} ;
    rname = '*' | [!-()+-<>-~] [!-~]* ;
    pos = [0-9]{1,9} ;
    mapq = [0-9]{1,3} ;

    action cigar_start_op_length { cigar_op_len_start = p - line.ptr; }
    action cigar_set_op_length {
        cigar_op_len = to!uint(line[cigar_op_len_start .. p - line.ptr]);
    }

    action cigar_set_op_chr { cigar_op_chr = fc; }
    action cigar_put_operation { 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }

    cigar = '*' | (digit+ > cigar_start_op_length % cigar_set_op_length
                  [MIDNSHPX=] > cigar_set_op_chr % cigar_put_operation)+ ;

    rnext = '*' | '=' | [!-()+-<>-~][!-~]* ;
    pnext = digit{1,9} ;
    tlen = sign? digit{1,9} ;

    action sequence_start { sequence_beg = p - line.ptr; }
    action sequence_end { sequence_end = p - line.ptr; }

    seq = '*' | ([A-Za-z=.]+ > sequence_start % sequence_end) ;
    qual = [!-~]+ ;

    mandatoryfields = qname '\t' 
                       flag '\t' 
                      rname '\t' 
                        pos '\t' 
                       mapq '\t'
                      cigar '\t'
                      rnext '\t'
                      pnext '\t'
                       tlen '\t'
                        seq '\t'
                       qual ;

    charvalue =  [!-~] ;
    integervalue = sign? digit+ ;
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

Alignment parseAlignmentLine(string line) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = null;
    int cs;

    size_t read_name_beg = 0;
    size_t read_name_end = 0;

    size_t sequence_beg = 0;
    size_t sequence_end = 0;

    uint cigar_op_len;
    char cigar_op_chr;

    size_t cigar_op_len_start;
    
    auto cigar = appender!(CigarOperation[])();

    %%write init;
    %%write exec;

    auto read = Alignment(line[read_name_beg .. read_name_end],
                          line[sequence_beg .. sequence_end],
                          cigar.data);

    return read;
}

unittest {
    import std.algorithm;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35 MD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";

	auto alignment = parseAlignmentLine(line);
    assert(alignment.read_name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
}
