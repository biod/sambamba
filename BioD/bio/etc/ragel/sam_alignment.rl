/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

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

    ### 1. STORE READ NAME ###
    action qname_start { read_name_beg = p - line.ptr; }
    action qname_end { read_name = line[read_name_beg .. p - line.ptr]; }
    action handle_invalid_qname { fhold; fgoto recover_from_invalid_qname; }
    recover_from_invalid_qname := invalid_field '\t' @{ fhold; fgoto flag_parsing; } ;

    qname =  '*' | (([!-?A-~]{1,255})** > qname_start % qname_end) ;

    ### 2. STORE FLAG ###
    action set_flag { flag = to!ushort(int_value); }
    flag = uint % set_flag;
    action handle_invalid_flag { fhold; fgoto recover_from_invalid_flag; }
    recover_from_invalid_flag := invalid_field '\t' @{ fhold; fgoto rname_parsing; } ;

    ### 3. STORE RNAME ###
    action rname_start { rname_beg = p - line.ptr; }
    action rname_end {
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]);
    }

    action handle_invalid_rname { fhold; fgoto recover_from_invalid_rname; }
    recover_from_invalid_rname := invalid_field '\t' @{ fhold; fgoto pos_parsing; } ;

    rname = '*' | (([!-()+-<>-~] [!-~]*) > rname_start % rname_end);

    ### 4. STORE POS ###
    action set_pos { end_pos = pos = to!uint(int_value); }
    pos = uint % set_pos;
    action handle_invalid_pos { fhold; fgoto recover_from_invalid_pos; }
    recover_from_invalid_pos := invalid_field '\t' @{ fhold; fgoto mapq_parsing; } ;

    ### 5. STORE MAPPING QUALITY ###
    action set_mapping_quality { mapping_quality = to!ubyte(int_value); }
    mapq = uint % set_mapping_quality;
    action handle_invalid_mapq { fhold; fgoto recover_from_invalid_mapq; }
    recover_from_invalid_mapq := invalid_field '\t' @{ fhold; fgoto cigar_parsing; } ;

    ### 6. INITIALIZE OUTPUT BUFFER ###
    action init_output_buffer {
        buffer.capacity = 32 + read_name.length + 1;
        buffer.putUnsafe!int(ref_id);
        buffer.putUnsafe!int(pos - 1);

        enforce(read_name.length + 1 <= 255, "Read name " ~ read_name ~ " is too long!");

        // bin will be set later
        auto bin_mq_nl = ((cast(uint)mapping_quality) << 8) | (read_name.length + 1);
        buffer.putUnsafe(cast(uint)bin_mq_nl);

        // number of CIGAR operations will be set later
        buffer.putUnsafe!uint(flag << 16);

        buffer.putUnsafe!int(0);
        buffer.putUnsafe!int(-1); // mate ref. id
        buffer.putUnsafe!int(-1); // mate pos
        buffer.putUnsafe!int(0);  // tlen

        buffer.putUnsafe(cast(ubyte[])read_name);
        buffer.putUnsafe!ubyte(0);

        rollback_size = buffer.length;
    }

    ### 7. STORE CIGAR OPERATIONS ###
    action cigar_set_op_length { cigar_op_len = to!uint(int_value); }
    action cigar_set_op_chr { cigar_op_chr = fc; }
    action cigar_put_operation {
        auto op = CigarOperation(cigar_op_len, cigar_op_chr);
        if (op.is_reference_consuming)
            end_pos += op.length;
        buffer.put!CigarOperation(op);
        {
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) + 1;
        }
    }

    action handle_invalid_cigar {
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) & 0xFFFF0000;
        buffer.shrink(rollback_size);
        end_pos = pos + 1;
        fhold; fgoto recover_from_invalid_cigar;
    }
    recover_from_invalid_cigar := invalid_field '\t' @{ fhold; fgoto rnext_parsing; } ;

    cigar = '*' | (uint % cigar_set_op_length
                  [MIDNSHPX=] > cigar_set_op_chr % cigar_put_operation)+ ;

    ### 8. SET BIN ###
    action set_bin {
        if (end_pos == pos)
            ++end_pos;
        {
        auto bin = reg2bin(pos - 1, end_pos - 1); // 0-based [) interval
        auto ptr = cast(uint*)(buffer.data.ptr + 2 * uint.sizeof);
        *ptr = (*ptr) | ((cast(uint)bin) << 16);
        }
    }

    ### 9. SET MATE REF. ID ###
    action set_same_mate_ref_id {
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = ref_id;
        }
    }

    action rnext_start { rnext_beg = p - line.ptr; }
    action rnext_end {
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
        }
    }
    action handle_invalid_rnext { fhold; fgoto recover_from_invalid_rnext; }
    recover_from_invalid_rnext := invalid_field '\t' @{ fhold; fgoto pnext_parsing; } ;

    rnext = '*' | ('=' % set_same_mate_ref_id) |
                  (([!-()+-<>-~][!-~]*) > rnext_start % rnext_end) ;

    ### 10. SET MATE POSITION ###
    action set_mate_pos {
        {
        auto ptr = cast(int*)(buffer.data.ptr + 6 * int.sizeof);
        *ptr = to!int(int_value) - 1;
        }
    }
    action handle_invalid_pnext { fhold; fgoto recover_from_invalid_pnext; }
    recover_from_invalid_pnext := invalid_field '\t' @{ fhold; fgoto tlen_parsing; } ;

    pnext = uint % set_mate_pos;

    ### 11. SET TEMPLATE LENGTH ###
    action set_template_length {
        {
        auto ptr = cast(int*)(buffer.data.ptr + 7 * int.sizeof);
        *ptr = to!int(int_value);
        }
    }
    action handle_invalid_tlen { fhold; fgoto recover_from_invalid_tlen; }
    recover_from_invalid_tlen := invalid_field '\t' @{ fhold; fgoto seq_parsing; } ;

    tlen = int % set_template_length;

    ### 12. SET SEQUENCE ###
    action sequence_start { sequence_beg = p - line.ptr; }
    action sequence_end {
        auto data = cast(ubyte[])line[sequence_beg .. p - line.ptr];
        l_seq = cast(int)data.length;
        auto raw_len = (l_seq + 1) / 2;

        // reserve space for base qualities, too
        buffer.capacity = buffer.length + raw_len + l_seq;

        for (size_t i = 0; i < raw_len; ++i) {
            auto b = cast(ubyte)(Base(data[2 * i]).internal_code << 4);
            if (2 * i + 1 < l_seq)
                b |= cast(ubyte)(Base(data[2 * i + 1]).internal_code);
            buffer.putUnsafe!ubyte(b);
        }

        // set l_seq
        {
        auto ptr = cast(int*)(buffer.data.ptr + 4 * int.sizeof);
        *ptr = l_seq;
        }

        rollback_size = buffer.length;
    }
    action handle_invalid_seq {
        rollback_size = buffer.length;
        fhold; fgoto recover_from_invalid_seq;
    }
    recover_from_invalid_seq := invalid_field '\t' @{ fhold; fgoto qual_parsing; } ;

    action set_rollback_size {
        rollback_size = buffer.length;
    }

    seq = ('*' % set_rollback_size) | ([A-Za-z=.]+ > sequence_start % sequence_end) ;

    ### 13. SET BASE QUALITIES ###
    action convert_next_character_to_prob {
        ++quals_length;
        quals_last_char = fc;
        buffer.putUnsafe!ubyte(cast(ubyte)(fc - 33));
    }

    action qual_end {
        // '*' may correspond either to a one-base long sequence
        // or to absence of information
        if (quals_length == 1 && quals_last_char == '*' && l_seq == 0)
            buffer.shrink(rollback_size);
    }

    action handle_invalid_qual {
        buffer.shrink(rollback_size);
        for (size_t i = 0; i < l_seq; ++i)
            buffer.putUnsafe!ubyte(0xFF);
        rollback_size = buffer.length;
        fhold; fgoto recover_from_invalid_qual;
    }
    # FIXME
    recover_from_invalid_qual := invalid_field '\t' @{ fhold; fgoto tag_parsing; } ;

    action check_qual_length {
        if (buffer.length - rollback_size != l_seq) {
            buffer.shrink(rollback_size);
            for (size_t i = 0; i < l_seq; ++i)
                buffer.putUnsafe!ubyte(0xFF);
        }
        rollback_size = buffer.length;
    }
    qual = [!-~]+ $ convert_next_character_to_prob % qual_end ;

    ###### PARSE MANDATORY FIELDS #######
    mandatoryfields = qname_parsing: (qname $!handle_invalid_qname)
                      flag_parsing: '\t' (flag $!handle_invalid_flag)
                      rname_parsing: '\t' (rname $!handle_invalid_rname)
                      pos_parsing: '\t' (pos $!handle_invalid_pos)
                      mapq_parsing: '\t' (mapq $!handle_invalid_mapq)
                      cigar_parsing: '\t' > init_output_buffer (cigar $!handle_invalid_cigar)
                      rnext_parsing: '\t' > set_bin (rnext $!handle_invalid_rnext)
                      pnext_parsing: '\t' (pnext $!handle_invalid_pnext)
                      tlen_parsing: '\t' (tlen $!handle_invalid_tlen)
                      seq_parsing: '\t' (seq $!handle_invalid_seq)
                      qual_parsing: '\t' (qual % check_qual_length $!handle_invalid_qual) ;

    ############ TAG PARSING ######

    action set_charvalue {
        buffer.capacity = buffer.length + 4;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('A');
        buffer.putUnsafe!char(fc);
    }

    action set_integervalue {
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        if (int_value < 0) {
            if (int_value >= byte.min) {
                buffer.putUnsafe!char('c');
                buffer.putUnsafe(cast(byte)int_value);
            } else if (int_value >= short.min) {
                buffer.putUnsafe!char('s');
                buffer.putUnsafe(cast(short)int_value);
            } else if (int_value >= int.min) {
                buffer.putUnsafe!char('i');
                buffer.putUnsafe(cast(int)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        } else {
            if (int_value <= ubyte.max) {
                buffer.putUnsafe!char('C');
                buffer.putUnsafe(cast(ubyte)int_value);
            } else if (int_value <= ushort.max) {
                buffer.putUnsafe!char('S');
                buffer.putUnsafe(cast(ushort)int_value);
            } else if (int_value <= uint.max) {
                buffer.putUnsafe!char('I');
                buffer.putUnsafe(cast(uint)int_value);
            } else {
                throw new Exception("integer out of range");
            }
        }
    }

    action start_tagvalue { tagvalue_beg = p - line.ptr; }

    action set_floatvalue {
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('f');
        buffer.putUnsafe!float(float_value);
    }

    action set_stringvalue {
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('Z');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }

    action set_hexstringvalue {
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('H');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }

    charvalue =  [!-~] > set_charvalue ;
    integervalue = int % set_integervalue;
    floatvalue = float % set_floatvalue ;

    action start_arrayvalue {
        arraytype = fc;
        buffer.capacity = buffer.length + 8;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('B');
        buffer.putUnsafe!char(arraytype);
        buffer.putUnsafe!uint(0);
        tag_array_length_offset = buffer.length - uint.sizeof;
    }

    action put_integer_to_array {
        // here, we assume that compiler is smart enough to move switch out of loop.
        switch (arraytype) {
            case 'c': buffer.put(to!byte(int_value)); break;
            case 'C': buffer.put(to!ubyte(int_value)); break;
            case 's': buffer.put(to!short(int_value)); break;
            case 'S': buffer.put(to!ushort(int_value)); break;
            case 'i': buffer.put(to!int(int_value)); break;
            case 'I': buffer.put(to!uint(int_value)); break;
            default: assert(0);
        }
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }

    action put_float_to_array {
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }

    stringvalue = [ !-~]+ > start_tagvalue % set_stringvalue ;
    hexstringvalue = xdigit+ > start_tagvalue % set_hexstringvalue ;
    integerarrayvalue = [cCsSiI] > start_arrayvalue (',' int % put_integer_to_array)+ ;
    floatarrayvalue = [f] > start_arrayvalue (',' float % put_float_to_array)+ ;
    arrayvalue = integerarrayvalue | floatarrayvalue ;

    tagvalue = ("A:" charvalue) |
               ("i:" integervalue) |
               ("f:" floatvalue) |
               ("Z:" stringvalue) |
               ("H:" hexstringvalue) |
               ("B:" arrayvalue) ;

    action tag_key_start { tag_key_beg = p - line.ptr; }
    action tag_key_end   { tag_key = cast(ubyte[])(line[tag_key_beg .. p - line.ptr]); }

    action handle_invalid_tag {
        buffer.shrink(rollback_size);
        fhold; fgoto recover_from_invalid_tag;
    }
    # FIXME: what if the tag is last?
    recover_from_invalid_tag := invalid_field '\t' @{ fhold; fgoto tag_parsing; } ;

    action update_rollback_size { rollback_size = buffer.length; }
    tag = (alpha alnum) > tag_key_start % tag_key_end ;
    optionalfield = tag ':' tagvalue % update_rollback_size $!handle_invalid_tag ;
    optionalfields = optionalfield ('\t' optionalfield)* ;

    alignment := field_parsing: mandatoryfields
                 tag_parsing: ('\t' optionalfields)? ;

    write data;
}%%

import bio.std.hts.sam.header;
import bio.std.hts.bam.cigar;
import bio.std.hts.bam.read;
import bio.std.hts.bam.bai.bin;
import bio.core.utils.outbuffer;
import bio.core.base;
import std.conv;
import std.array;
import std.exception;

BamRead parseAlignmentLine(string line, SamHeader header, OutBuffer buffer=null) {
    char* p = cast(char*)line.ptr;
    char* pe = p + line.length;
    char* eof = pe;
    int cs;

    if (buffer is null)
        buffer = new OutBuffer(8192);
    else
        buffer.clear();

    size_t rollback_size; // needed in case of invalid data

    byte current_sign = 1;

    size_t read_name_beg; // position of beginning of QNAME

    size_t sequence_beg; // position of SEQ start
    int l_seq;           // sequence length

    uint cigar_op_len;   // length of CIGAR operation
    char cigar_op_chr;   // CIGAR operation

    size_t quals_length;  // number of QUAL characters
    char quals_last_char; // needed in order to handle '*' correctly

    size_t cigar_op_len_start; // position of start of CIGAR operation

    long int_value;                      // for storing temporary integers
    float float_value;                   // for storing temporary floats
    size_t float_beg;                    // position of start of current float
    char arraytype;                      // type of last array tag value
    size_t tag_array_length_offset;      // where the length is stored in the buffer

    string read_name;
    ushort flag;
    int pos = -1;
    int end_pos; // for bin calculation
    int mate_pos = -1;
    ubyte mapping_quality = 255;
    int template_length = 0;

    size_t tag_key_beg, tagvalue_beg;
    ubyte[] tag_key;
    size_t rname_beg, rnext_beg;

    int ref_id = -1;

    %%write init;
    %%write exec;

    BamRead read;
    read.raw_data = buffer.data[];
    return read;
}

unittest {
    import std.algorithm;
    import std.math;

    auto line = "ERR016155.15021091\t185\t20\t60033\t25\t66S35M\t=\t60033\t0\tAGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC\t#####################################################################################################\tX0:i:1\tX1:i:0\tXC:i:35\tMD:Z:17A8A8\tRG:Z:ERR016155\tAM:i:0\tNM:i:2\tSM:i:25\tXT:A:U\tBQ:Z:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\tY0:B:c,1,2,3\tY1:B:f,13.263,-3.1415,52.63461";

    auto header = new SamHeader("@SQ\tSN:20\tLN:1234567");
    auto alignment = parseAlignmentLine(line, header);
    assert(alignment.name == "ERR016155.15021091");
    assert(equal(alignment.sequence(), "AGAAAAAACTGGAAGTTAATAGAGTGGTGACTCAGATCCAGTGGTGGAAGGGTAAGGGATCTTGGAACCCTATAGAGTTGCTGTGTGCCAGGGCCAGATCC"));
    assert(alignment.cigarString() == "66S35M");
    assert(alignment.flag == 185);
    assert(alignment.position == 60032);
    assert(alignment.mapping_quality == 25);
    assert(alignment.mate_position == 60032);
    assert(alignment.ref_id == 0);
    assert(alignment.mate_ref_id == 0);
    assert(to!ubyte(alignment["AM"]) == 0);
    assert(to!ubyte(alignment["SM"]) == 25);
    assert(to!string(alignment["MD"]) == "17A8A8");
    assert(equal(to!(byte[])(alignment["Y0"]), [1, 2, 3]));
    assert(equal!approxEqual(to!(float[])(alignment["Y1"]), [13.263, -3.1415, 52.63461]));
    assert(to!char(alignment["XT"]) == 'U');

    import bio.std.hts.bam.reference;

    auto info = ReferenceSequenceInfo("20", 1234567);

    auto invalid_cigar_string = "1\t100\t20\t50000\t30\tMZABC\t=\t50000\t0\tACGT\t####";
    alignment = parseAlignmentLine(invalid_cigar_string, header);
    assert(equal(alignment.sequence(), "ACGT"));

    auto invalid_tag_and_qual = "2\t100\t20\t5\t40\t27M30X5D\t=\t3\t10\tACT\t !\n\tX1:i:7\tX3:i:zzz\tX4:i:5";
    alignment = parseAlignmentLine(invalid_tag_and_qual, header);
    assert(alignment.base_qualities == [255, 255, 255]); // i.e. invalid
    assert(to!ubyte(alignment["X1"]) == 7);
    assert(alignment["X3"].is_nothing);
    assert(to!ubyte(alignment["X4"]) == 5);
}
