/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
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

    float = ((sign? ((digit* '.'? digit+ ([eE] sign? digit+)?) | "inf") ) | "nan")
                > mark_float_start % update_float_value ;

    invalid_field = [^\t]* ; # TODO: make class with onError methods and pass it to parseAlignmentLine

    action qname_start { read_name_beg = p - line.ptr; }
    action qname_end { read_name_end = p - line.ptr; }

    qname =  '*' | (([!-?A-~]{1,255})** > qname_start % qname_end) ;

    action set_flag { flag = to!ushort(int_value); }
    action set_pos { pos = to!uint(int_value); }
    action set_mapping_quality { mapping_quality = to!ubyte(int_value); }

    flag = uint % set_flag;

    action rname_start { rname_beg = p - line.ptr; }
    action rname_end {
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]); 
    }

    rname = '*' | (([!-()+-<>-~] [!-~]*) > rname_start % rname_end);

    pos = uint % set_pos;
    mapq = uint % set_mapping_quality;

    action cigar_set_op_length { cigar_op_len = to!uint(int_value); }
    action cigar_set_op_chr { cigar_op_chr = fc; }
    action cigar_put_operation { 
        cigar.put(CigarOperation(cigar_op_len, cigar_op_chr)); 
    }

    cigar = '*' | (uint % cigar_set_op_length
                  [MIDNSHPX=] > cigar_set_op_chr % cigar_put_operation)+ ;

    action set_same_mate_ref_id {
        mate_ref_id = ref_id;
    }
   
    action rnext_start { rnext_beg = p - line.ptr; }
    action rnext_end {
        mate_ref_id = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
    }

    rnext = '*' | ('=' % set_same_mate_ref_id) | 
                  (([!-()+-<>-~][!-~]*) > rnext_start % rnext_end) ;

    action set_mate_pos { mate_pos = to!uint(int_value); }
    action set_template_length { template_length = to!int(int_value); }

    pnext = uint % set_mate_pos;
    tlen = int % set_template_length;

    action sequence_start { sequence_beg = p - line.ptr; }
    action sequence_end { sequence = line[sequence_beg .. p - line.ptr]; }

    seq = '*' | ([A-Za-z=.]+ > sequence_start % sequence_end) ;
    qual = [!-~]+ ;

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

    mandatoryfields = (qname | invalid_field) '\t'
                      (flag  | invalid_field) '\t'
                      (rname | invalid_field) '\t'
                      (pos   | invalid_field) '\t'
                      (mapq  | invalid_field) '\t'
                      (cigar | invalid_field) '\t'
                      (rnext | invalid_field) '\t'
                      (pnext | invalid_field) '\t'
                      (tlen  | invalid_field) '\t'
                      (seq   | invalid_field) '\t'
                      ((qual > allocate_quality_array
                             $ convert_next_character_to_prob) | invalid_field) ;

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
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
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
    action append_tag_value { builder.put(current_tag, current_tagvalue); }

    tag = (alpha alnum) > tag_key_start % tag_key_end ;
    optionalfield = (tag ':' tagvalue % append_tag_value) | invalid_field ;
    optionalfields = optionalfield ('\t' optionalfield)* ;

    alignment := mandatoryfields ('\t' optionalfields)? ;

    write data; 
}%%

import alignment;
import tagvalue;
import samheader;

import std.array;
import std.conv;
import std.typecons;
import std.outbuffer;
import std.c.stdlib;
import utils.tagstoragebuilder;

class AlignmentBuildStorage {
    Appender!(CigarOperation[]) cigar_appender;
    OutBuffer outbuffer;
    TagStorageBuilder tag_storage_builder;

    this() {
        cigar_appender = appender!(CigarOperation[])();
        outbuffer = new OutBuffer();
        tag_storage_builder = TagStorageBuilder.create();
    }

    void clear() {
        cigar_appender.clear();
        tag_storage_builder.clear();
        outbuffer.data.length = 0;
        outbuffer.offset = 0;
    }
}

Alignment parseAlignmentLine(string line, SamHeader header, 
                             AlignmentBuildStorage b=null) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    if (b is null) {
        b = new AlignmentBuildStorage();
    } else {
        b.clear();
    }

    byte current_sign = 1;

    size_t read_name_beg; // position of beginning of QNAME
    size_t read_name_end; // position past the end of QNAME

    size_t sequence_beg; // position of SEQ start
    string sequence;     // SEQ

    uint cigar_op_len;   // length of CIGAR operation
    char cigar_op_chr;   // CIGAR operation

    size_t cigar_op_len_start; // position of start of CIGAR operation
    
    auto cigar = b.cigar_appender;

    long int_value;                      // for storing temporary integers
    float float_value;                   // for storing temporary floats
    size_t float_beg;                    // position of start of current float
    auto outbuffer = b.outbuffer;        // used to build tag values which hold arrays
    char arraytype;                      // type of last array tag value

    ushort flag;
    uint pos;
    uint mate_pos;
    ubyte mapping_quality; 
    int template_length;
    ubyte* qual_ptr = null;
    size_t qual_index; 

    string current_tag;
    Value current_tagvalue;

    size_t tag_key_beg, tagvalue_beg;
    size_t rname_beg, rnext_beg;

    int ref_id = -1;
    int mate_ref_id = -1;
    
    auto builder = b.tag_storage_builder;

    %%write init;
    %%write exec;

    auto read = Alignment(line[read_name_beg .. read_name_end], 
                          sequence,
                          cigar.data,
                          builder.data);

    if (qual_ptr !is null && qual_index == sequence.length) {
        read.phred_base_quality = qual_ptr[0 .. sequence.length];
    }

    read.flag = flag;
    read.mapping_quality = mapping_quality;
    read.position = pos - 1; // we use 0-based coordinates, not 1-based
    read.template_length = template_length;
    read.next_pos = mate_pos - 1; // also 0-based
    read.ref_id = ref_id;
    read.next_ref_id = mate_ref_id;

    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto header = new SamHeader("@SQ\tSN:20\tLN:1234567");
    auto alignment = parseAlignmentLine(line, header);
    assert(alignment.read_name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60032);
    assert(alignment.mapping_quality == 25);
    assert(alignment.next_pos == 60032);
    assert(alignment.ref_id == 0);
    assert(alignment.next_ref_id == 0);
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
    assert(to!char(alignment["XT"]) == 'U');

    import std.stdio;
    import sam.serialize;
    import reference;

    ReferenceSequenceInfo info;
    info.name = "20";
    info.length = 1234567;

    auto invalid_cigar_string = "1\t100\t20\t50000\t30\tMZABC\t=\t50000\t0\tACGT\t####";
    alignment = parseAlignmentLine(invalid_cigar_string, header);
    assert(equal(alignment.sequence(), "ACGT"));

    auto invalid_tag_and_qual = "2\t100\t20\t5\t40\t27M30X5D\t=\t3\t10\tACT\t !\n\tX1:i:7\tX3:i:zzz\tX4:i:5";
    alignment = parseAlignmentLine(invalid_tag_and_qual, header);
    assert(alignment.phred_base_quality == [255, 255, 255]); // i.e. invalid
    assert(to!ubyte(alignment["X1"]) == 7);
    assert(alignment["X3"].is_nothing);
    assert(to!ubyte(alignment["X4"]) == 5);

}
