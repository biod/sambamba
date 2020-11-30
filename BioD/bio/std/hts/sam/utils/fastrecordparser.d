module bio.std.hts.sam.utils.fastrecordparser;

#line 1 "sam_alignment.rl"
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

#line 28 "sam_alignment.d"
static const int sam_alignment_start = 1;
static const int sam_alignment_first_final = 191;
static const int sam_alignment_error = 0;

static const int sam_alignment_en_recover_from_invalid_qname = 169;
static const int sam_alignment_en_recover_from_invalid_flag = 170;
static const int sam_alignment_en_recover_from_invalid_rname = 171;
static const int sam_alignment_en_recover_from_invalid_pos = 172;
static const int sam_alignment_en_recover_from_invalid_mapq = 173;
static const int sam_alignment_en_recover_from_invalid_cigar = 174;
static const int sam_alignment_en_recover_from_invalid_rnext = 175;
static const int sam_alignment_en_recover_from_invalid_pnext = 176;
static const int sam_alignment_en_recover_from_invalid_tlen = 177;
static const int sam_alignment_en_recover_from_invalid_seq = 178;
static const int sam_alignment_en_recover_from_invalid_qual = 179;
static const int sam_alignment_en_recover_from_invalid_tag = 180;
static const int sam_alignment_en_alignment = 1;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_flag_parsing = 181;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_rname_parsing = 182;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_pos_parsing = 183;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_mapq_parsing = 184;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_cigar_parsing = 185;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_rnext_parsing = 186;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_pnext_parsing = 187;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_tlen_parsing = 188;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_seq_parsing = 189;
static const int sam_alignment_en_alignment_field_parsing_mandatoryfields_qual_parsing = 190;
static const int sam_alignment_en_alignment_tag_parsing = 251;


#line 419 "sam_alignment.rl"


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

    
#line 121 "sam_alignment.d"
	{
	cs = sam_alignment_start;
	}

#line 480 "sam_alignment.rl"
    
#line 128 "sam_alignment.d"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
goto case; case 1:
	if ( (*p) == 9u )
		goto tr1;
	if ( (*p) > 63u ) {
		if ( 65u <= (*p) && (*p) <= 126u )
			goto tr2;
	} else if ( (*p) >= 33u )
		goto tr2;
	goto tr0;
tr0:
#line 50 "sam_alignment.rl"
	{ p--; {if (true) goto st169;} }
	goto st0;
tr3:
#line 58 "sam_alignment.rl"
	{ p--; {if (true) goto st170;} }
	goto st0;
tr7:
#line 67 "sam_alignment.rl"
	{ p--; {if (true) goto st171;} }
	goto st0;
tr12:
#line 75 "sam_alignment.rl"
	{ p--; {if (true) goto st172;} }
	goto st0;
tr16:
#line 81 "sam_alignment.rl"
	{ p--; {if (true) goto st173;} }
	goto st0;
tr20:
#line 124 "sam_alignment.rl"
	{
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) & 0xFFFF0000;
        buffer.shrink(rollback_size);
        end_pos = pos + 1;
        p--; {if (true) goto st174;}
    }
	goto st0;
tr24:
#line 162 "sam_alignment.rl"
	{ p--; {if (true) goto st175;} }
	goto st0;
tr30:
#line 175 "sam_alignment.rl"
	{ p--; {if (true) goto st176;} }
	goto st0;
tr34:
#line 187 "sam_alignment.rl"
	{ p--; {if (true) goto st177;} }
	goto st0;
tr39:
#line 217 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
        p--; {if (true) goto st178;}
    }
	goto st0;
tr43:
#line 243 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        for (size_t i = 0; i < l_seq; ++i)
            buffer.putUnsafe!ubyte(0xFF);
        rollback_size = buffer.length;
        p--; {if (true) goto st179;}
    }
	goto st0;
tr49:
#line 403 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        p--; {if (true) goto st180;}
    }
	goto st0;
#line 209 "sam_alignment.d"
st0:
cs = 0;
	goto _out;
tr1:
#line 48 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
#line 49 "sam_alignment.rl"
	{ read_name = line[read_name_beg .. p - line.ptr]; }
	goto st2;
tr206:
#line 49 "sam_alignment.rl"
	{ read_name = line[read_name_beg .. p - line.ptr]; }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
goto case; case 2:
#line 227 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr4;
	goto tr3;
tr4:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
goto case; case 3:
#line 241 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr6;
	goto tr3;
tr5:
#line 56 "sam_alignment.rl"
	{ flag = to!ushort(int_value); }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
goto case; case 4:
#line 255 "sam_alignment.d"
	if ( (*p) == 42u )
		goto st150;
	if ( (*p) > 60u ) {
		if ( 62u <= (*p) && (*p) <= 126u )
			goto tr8;
	} else if ( (*p) >= 33u )
		goto tr8;
	goto tr7;
tr8:
#line 62 "sam_alignment.rl"
	{ rname_beg = p - line.ptr; }
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
goto case; case 5:
#line 272 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr10;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st5;
	goto tr7;
tr10:
#line 63 "sam_alignment.rl"
	{
        ref_id = header.getSequenceIndex(line[rname_beg .. p - line.ptr]);
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
goto case; case 6:
#line 288 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr13;
	goto tr12;
tr13:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
goto case; case 7:
#line 302 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr15;
	goto tr12;
tr14:
#line 73 "sam_alignment.rl"
	{ end_pos = pos = to!uint(int_value); }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
goto case; case 8:
#line 316 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr17;
	goto tr16;
tr17:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
goto case; case 9:
#line 330 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr19;
	goto tr16;
tr18:
#line 79 "sam_alignment.rl"
	{ mapping_quality = to!ubyte(int_value); }
#line 85 "sam_alignment.rl"
	{
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
	goto st10;
tr235:
#line 85 "sam_alignment.rl"
	{
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
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
goto case; case 10:
#line 396 "sam_alignment.d"
	if ( (*p) == 42u )
		goto st11;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr22;
	goto tr20;
st11:
	if ( ++p == pe )
		goto _test_eof11;
goto case; case 11:
	if ( (*p) == 9u )
		goto tr23;
	goto tr20;
tr23:
#line 137 "sam_alignment.rl"
	{
        if (end_pos == pos)
            ++end_pos;
        {
        auto bin = reg2bin(pos - 1, end_pos - 1); // 0-based [) interval
        auto ptr = cast(uint*)(buffer.data.ptr + 2 * uint.sizeof);
        *ptr = (*ptr) | ((cast(uint)bin) << 16);
        }
    }
	goto st12;
tr155:
#line 113 "sam_alignment.rl"
	{
        auto op = CigarOperation(cigar_op_len, cigar_op_chr);
        if (op.is_reference_consuming)
            end_pos += op.length;
        buffer.put!CigarOperation(op);
        {
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) + 1;
        }
    }
#line 137 "sam_alignment.rl"
	{
        if (end_pos == pos)
            ++end_pos;
        {
        auto bin = reg2bin(pos - 1, end_pos - 1); // 0-based [) interval
        auto ptr = cast(uint*)(buffer.data.ptr + 2 * uint.sizeof);
        *ptr = (*ptr) | ((cast(uint)bin) << 16);
        }
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
goto case; case 12:
#line 448 "sam_alignment.d"
	switch( (*p) ) {
		case 42u: goto st95;
		case 61u: goto st96;
		default: break;
	}
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr25;
	goto tr24;
tr25:
#line 155 "sam_alignment.rl"
	{ rnext_beg = p - line.ptr; }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
goto case; case 13:
#line 465 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr28;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st13;
	goto tr24;
tr28:
#line 156 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = header.getSequenceIndex(line[rnext_beg .. p - line.ptr]);
        }
    }
	goto st14;
tr136:
#line 148 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 5 * int.sizeof);
        *ptr = ref_id;
        }
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
goto case; case 14:
#line 493 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr31;
	goto tr30;
tr31:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
goto case; case 15:
#line 507 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr33;
	goto tr30;
tr32:
#line 169 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 6 * int.sizeof);
        *ptr = to!int(int_value) - 1;
        }
    }
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
goto case; case 16:
#line 526 "sam_alignment.d"
	switch( (*p) ) {
		case 43u: goto tr35;
		case 45u: goto tr35;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr36;
	goto tr34;
tr35:
#line 27 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
goto case; case 17:
#line 543 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr36;
	goto tr34;
tr36:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
goto case; case 18:
#line 557 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr38;
	goto tr34;
tr37:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 181 "sam_alignment.rl"
	{
        {
        auto ptr = cast(int*)(buffer.data.ptr + 7 * int.sizeof);
        *ptr = to!int(int_value);
        }
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
goto case; case 19:
#line 578 "sam_alignment.d"
	switch( (*p) ) {
		case 42u: goto st20;
		case 46u: goto tr41;
		case 61u: goto tr41;
		default: break;
	}
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto tr41;
	} else if ( (*p) >= 65u )
		goto tr41;
	goto tr39;
st20:
	if ( ++p == pe )
		goto _test_eof20;
goto case; case 20:
	if ( (*p) == 9u )
		goto tr42;
	goto tr39;
tr42:
#line 223 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
    }
	goto st21;
tr101:
#line 194 "sam_alignment.rl"
	{
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
	goto st21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
goto case; case 21:
#line 634 "sam_alignment.d"
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr44;
	goto tr43;
tr44:
#line 230 "sam_alignment.rl"
	{
        ++quals_length;
        quals_last_char = (*p);
        buffer.putUnsafe!ubyte(cast(ubyte)((*p) - 33));
    }
	goto st191;
st191:
	if ( ++p == pe )
		goto _test_eof191;
goto case; case 191:
#line 650 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr239;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr44;
	goto tr43;
tr239:
#line 236 "sam_alignment.rl"
	{
        // '*' may correspond either to a one-base long sequence
        // or to absence of information
        if (quals_length == 1 && quals_last_char == '*' && l_seq == 0)
            buffer.shrink(rollback_size);
    }
#line 253 "sam_alignment.rl"
	{
        if (buffer.length - rollback_size != l_seq) {
            buffer.shrink(rollback_size);
            for (size_t i = 0; i < l_seq; ++i)
                buffer.putUnsafe!ubyte(0xFF);
        }
        rollback_size = buffer.length;
    }
	goto st22;
tr240:
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr241:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 362 "sam_alignment.rl"
	{
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
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr260:
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 379 "sam_alignment.rl"
	{
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr263:
#line 337 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('H');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr265:
#line 326 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('Z');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr267:
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 319 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('f');
        buffer.putUnsafe!float(float_value);
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
tr269:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 285 "sam_alignment.rl"
	{
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
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
goto case; case 22:
#line 804 "sam_alignment.d"
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto tr45;
	} else if ( (*p) >= 65u )
		goto tr45;
	goto st0;
tr45:
#line 400 "sam_alignment.rl"
	{ tag_key_beg = p - line.ptr; }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
goto case; case 23:
#line 819 "sam_alignment.d"
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto st24;
	} else if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto st24;
	} else
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
goto case; case 24:
	if ( (*p) == 58u )
		goto tr48;
	goto st0;
tr48:
#line 401 "sam_alignment.rl"
	{ tag_key = cast(ubyte[])(line[tag_key_beg .. p - line.ptr]); }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
goto case; case 25:
#line 844 "sam_alignment.d"
	switch( (*p) ) {
		case 65u: goto st26;
		case 66u: goto st28;
		case 72u: goto st43;
		case 90u: goto st45;
		case 102u: goto st47;
		case 105u: goto st57;
		default: break;
	}
	goto tr49;
st26:
	if ( ++p == pe )
		goto _test_eof26;
goto case; case 26:
	if ( (*p) == 58u )
		goto st27;
	goto tr49;
st27:
	if ( ++p == pe )
		goto _test_eof27;
goto case; case 27:
	if ( 33u <= (*p) && (*p) <= 126u )
		goto tr57;
	goto tr49;
tr57:
#line 278 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 4;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('A');
        buffer.putUnsafe!char((*p));
    }
	goto st192;
st192:
	if ( ++p == pe )
		goto _test_eof192;
goto case; case 192:
#line 882 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr240;
	goto tr49;
st28:
	if ( ++p == pe )
		goto _test_eof28;
goto case; case 28:
	if ( (*p) == 58u )
		goto st29;
	goto tr49;
st29:
	if ( ++p == pe )
		goto _test_eof29;
goto case; case 29:
	switch( (*p) ) {
		case 67u: goto tr59;
		case 73u: goto tr59;
		case 83u: goto tr59;
		case 99u: goto tr59;
		case 102u: goto tr60;
		case 105u: goto tr59;
		case 115u: goto tr59;
		default: break;
	}
	goto tr49;
tr59:
#line 352 "sam_alignment.rl"
	{
        arraytype = (*p);
        buffer.capacity = buffer.length + 8;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('B');
        buffer.putUnsafe!char(arraytype);
        buffer.putUnsafe!uint(0);
        tag_array_length_offset = buffer.length - uint.sizeof;
    }
	goto st30;
st30:
	if ( ++p == pe )
		goto _test_eof30;
goto case; case 30:
#line 924 "sam_alignment.d"
	if ( (*p) == 44u )
		goto st31;
	goto tr49;
tr242:
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 362 "sam_alignment.rl"
	{
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
	goto st31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
goto case; case 31:
#line 953 "sam_alignment.d"
	switch( (*p) ) {
		case 43u: goto tr62;
		case 45u: goto tr62;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr63;
	goto tr49;
tr62:
#line 27 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st32;
st32:
	if ( ++p == pe )
		goto _test_eof32;
goto case; case 32:
#line 970 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr63;
	goto tr49;
tr63:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st193;
st193:
	if ( ++p == pe )
		goto _test_eof193;
goto case; case 193:
#line 984 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr243;
	goto tr49;
tr243:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st194;
st194:
	if ( ++p == pe )
		goto _test_eof194;
goto case; case 194:
#line 1001 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr244;
	goto tr49;
tr244:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st195;
st195:
	if ( ++p == pe )
		goto _test_eof195;
goto case; case 195:
#line 1018 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr245;
	goto tr49;
tr245:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st196;
st196:
	if ( ++p == pe )
		goto _test_eof196;
goto case; case 196:
#line 1035 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr246;
	goto tr49;
tr246:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st197;
st197:
	if ( ++p == pe )
		goto _test_eof197;
goto case; case 197:
#line 1052 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr247;
	goto tr49;
tr247:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st198;
st198:
	if ( ++p == pe )
		goto _test_eof198;
goto case; case 198:
#line 1069 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr248;
	goto tr49;
tr248:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st199;
st199:
	if ( ++p == pe )
		goto _test_eof199;
goto case; case 199:
#line 1086 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr249;
	goto tr49;
tr249:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st200;
st200:
	if ( ++p == pe )
		goto _test_eof200;
goto case; case 200:
#line 1103 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr250;
	goto tr49;
tr250:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st201;
st201:
	if ( ++p == pe )
		goto _test_eof201;
goto case; case 201:
#line 1120 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr251;
	goto tr49;
tr251:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st202;
st202:
	if ( ++p == pe )
		goto _test_eof202;
goto case; case 202:
#line 1137 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr252;
	goto tr49;
tr252:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st203;
st203:
	if ( ++p == pe )
		goto _test_eof203;
goto case; case 203:
#line 1154 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr253;
	goto tr49;
tr253:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st204;
st204:
	if ( ++p == pe )
		goto _test_eof204;
goto case; case 204:
#line 1171 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr254;
	goto tr49;
tr254:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st205;
st205:
	if ( ++p == pe )
		goto _test_eof205;
goto case; case 205:
#line 1188 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr255;
	goto tr49;
tr255:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st206;
st206:
	if ( ++p == pe )
		goto _test_eof206;
goto case; case 206:
#line 1205 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr256;
	goto tr49;
tr256:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st207;
st207:
	if ( ++p == pe )
		goto _test_eof207;
goto case; case 207:
#line 1222 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr257;
	goto tr49;
tr257:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st208;
st208:
	if ( ++p == pe )
		goto _test_eof208;
goto case; case 208:
#line 1239 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr258;
	goto tr49;
tr258:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st209;
st209:
	if ( ++p == pe )
		goto _test_eof209;
goto case; case 209:
#line 1256 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr259;
	goto tr49;
tr259:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st210;
st210:
	if ( ++p == pe )
		goto _test_eof210;
goto case; case 210:
#line 1273 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr241;
		case 44u: goto tr242;
		default: break;
	}
	goto tr49;
tr60:
#line 352 "sam_alignment.rl"
	{
        arraytype = (*p);
        buffer.capacity = buffer.length + 8;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('B');
        buffer.putUnsafe!char(arraytype);
        buffer.putUnsafe!uint(0);
        tag_array_length_offset = buffer.length - uint.sizeof;
    }
	goto st33;
st33:
	if ( ++p == pe )
		goto _test_eof33;
goto case; case 33:
#line 1296 "sam_alignment.d"
	if ( (*p) == 44u )
		goto st34;
	goto tr49;
tr261:
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 379 "sam_alignment.rl"
	{
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
	goto st34;
st34:
	if ( ++p == pe )
		goto _test_eof34;
goto case; case 34:
#line 1318 "sam_alignment.d"
	switch( (*p) ) {
		case 43u: goto tr65;
		case 45u: goto tr65;
		case 46u: goto tr66;
		case 105u: goto tr68;
		case 110u: goto tr69;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr67;
	goto tr49;
tr65:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st35;
st35:
	if ( ++p == pe )
		goto _test_eof35;
goto case; case 35:
#line 1338 "sam_alignment.d"
	switch( (*p) ) {
		case 46u: goto st36;
		case 105u: goto st39;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st213;
	goto tr49;
tr66:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
goto case; case 36:
#line 1355 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st211;
	goto tr49;
st211:
	if ( ++p == pe )
		goto _test_eof211;
goto case; case 211:
	switch( (*p) ) {
		case 9u: goto tr260;
		case 44u: goto tr261;
		case 69u: goto st37;
		case 101u: goto st37;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st211;
	goto tr49;
st37:
	if ( ++p == pe )
		goto _test_eof37;
goto case; case 37:
	switch( (*p) ) {
		case 43u: goto st38;
		case 45u: goto st38;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st212;
	goto tr49;
st38:
	if ( ++p == pe )
		goto _test_eof38;
goto case; case 38:
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st212;
	goto tr49;
st212:
	if ( ++p == pe )
		goto _test_eof212;
goto case; case 212:
	switch( (*p) ) {
		case 9u: goto tr260;
		case 44u: goto tr261;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st212;
	goto tr49;
tr67:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st213;
st213:
	if ( ++p == pe )
		goto _test_eof213;
goto case; case 213:
#line 1412 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr260;
		case 44u: goto tr261;
		case 46u: goto st36;
		case 69u: goto st37;
		case 101u: goto st37;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st213;
	goto tr49;
tr68:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st39;
st39:
	if ( ++p == pe )
		goto _test_eof39;
goto case; case 39:
#line 1432 "sam_alignment.d"
	if ( (*p) == 110u )
		goto st40;
	goto tr49;
st40:
	if ( ++p == pe )
		goto _test_eof40;
goto case; case 40:
	if ( (*p) == 102u )
		goto st214;
	goto tr49;
st214:
	if ( ++p == pe )
		goto _test_eof214;
goto case; case 214:
	switch( (*p) ) {
		case 9u: goto tr260;
		case 44u: goto tr261;
		default: break;
	}
	goto tr49;
tr69:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
goto case; case 41:
#line 1461 "sam_alignment.d"
	if ( (*p) == 97u )
		goto st42;
	goto tr49;
st42:
	if ( ++p == pe )
		goto _test_eof42;
goto case; case 42:
	if ( (*p) == 110u )
		goto st214;
	goto tr49;
st43:
	if ( ++p == pe )
		goto _test_eof43;
goto case; case 43:
	if ( (*p) == 58u )
		goto st44;
	goto tr49;
st44:
	if ( ++p == pe )
		goto _test_eof44;
goto case; case 44:
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr80;
	} else if ( (*p) > 70u ) {
		if ( 97u <= (*p) && (*p) <= 102u )
			goto tr80;
	} else
		goto tr80;
	goto tr49;
tr80:
#line 317 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	goto st215;
st215:
	if ( ++p == pe )
		goto _test_eof215;
goto case; case 215:
#line 1500 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr263;
	if ( (*p) < 65u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto st215;
	} else if ( (*p) > 70u ) {
		if ( 97u <= (*p) && (*p) <= 102u )
			goto st215;
	} else
		goto st215;
	goto tr49;
st45:
	if ( ++p == pe )
		goto _test_eof45;
goto case; case 45:
	if ( (*p) == 58u )
		goto st46;
	goto tr49;
st46:
	if ( ++p == pe )
		goto _test_eof46;
goto case; case 46:
	if ( 32u <= (*p) && (*p) <= 126u )
		goto tr82;
	goto tr49;
tr82:
#line 317 "sam_alignment.rl"
	{ tagvalue_beg = p - line.ptr; }
	goto st216;
st216:
	if ( ++p == pe )
		goto _test_eof216;
goto case; case 216:
#line 1534 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr265;
	if ( 32u <= (*p) && (*p) <= 126u )
		goto st216;
	goto tr49;
st47:
	if ( ++p == pe )
		goto _test_eof47;
goto case; case 47:
	if ( (*p) == 58u )
		goto st48;
	goto tr49;
st48:
	if ( ++p == pe )
		goto _test_eof48;
goto case; case 48:
	switch( (*p) ) {
		case 43u: goto tr84;
		case 45u: goto tr84;
		case 46u: goto tr85;
		case 105u: goto tr87;
		case 110u: goto tr88;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr86;
	goto tr49;
tr84:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
goto case; case 49:
#line 1570 "sam_alignment.d"
	switch( (*p) ) {
		case 46u: goto st50;
		case 105u: goto st53;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st219;
	goto tr49;
tr85:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
goto case; case 50:
#line 1587 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st217;
	goto tr49;
st217:
	if ( ++p == pe )
		goto _test_eof217;
goto case; case 217:
	switch( (*p) ) {
		case 9u: goto tr267;
		case 69u: goto st51;
		case 101u: goto st51;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st217;
	goto tr49;
st51:
	if ( ++p == pe )
		goto _test_eof51;
goto case; case 51:
	switch( (*p) ) {
		case 43u: goto st52;
		case 45u: goto st52;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st218;
	goto tr49;
st52:
	if ( ++p == pe )
		goto _test_eof52;
goto case; case 52:
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st218;
	goto tr49;
st218:
	if ( ++p == pe )
		goto _test_eof218;
goto case; case 218:
	if ( (*p) == 9u )
		goto tr267;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st218;
	goto tr49;
tr86:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st219;
st219:
	if ( ++p == pe )
		goto _test_eof219;
goto case; case 219:
#line 1640 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr267;
		case 46u: goto st50;
		case 69u: goto st51;
		case 101u: goto st51;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto st219;
	goto tr49;
tr87:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
goto case; case 53:
#line 1659 "sam_alignment.d"
	if ( (*p) == 110u )
		goto st54;
	goto tr49;
st54:
	if ( ++p == pe )
		goto _test_eof54;
goto case; case 54:
	if ( (*p) == 102u )
		goto st220;
	goto tr49;
st220:
	if ( ++p == pe )
		goto _test_eof220;
goto case; case 220:
	if ( (*p) == 9u )
		goto tr267;
	goto tr49;
tr88:
#line 37 "sam_alignment.rl"
	{ float_beg = p - line.ptr; }
	goto st55;
st55:
	if ( ++p == pe )
		goto _test_eof55;
goto case; case 55:
#line 1685 "sam_alignment.d"
	if ( (*p) == 97u )
		goto st56;
	goto tr49;
st56:
	if ( ++p == pe )
		goto _test_eof56;
goto case; case 56:
	if ( (*p) == 110u )
		goto st220;
	goto tr49;
st57:
	if ( ++p == pe )
		goto _test_eof57;
goto case; case 57:
	if ( (*p) == 58u )
		goto st58;
	goto tr49;
st58:
	if ( ++p == pe )
		goto _test_eof58;
goto case; case 58:
	switch( (*p) ) {
		case 43u: goto tr99;
		case 45u: goto tr99;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr100;
	goto tr49;
tr99:
#line 27 "sam_alignment.rl"
	{ current_sign = (*p) == '-' ? -1 : 1; }
	goto st59;
st59:
	if ( ++p == pe )
		goto _test_eof59;
goto case; case 59:
#line 1723 "sam_alignment.d"
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr100;
	goto tr49;
tr100:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st221;
st221:
	if ( ++p == pe )
		goto _test_eof221;
goto case; case 221:
#line 1737 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr270;
	goto tr49;
tr270:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st222;
st222:
	if ( ++p == pe )
		goto _test_eof222;
goto case; case 222:
#line 1751 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr271;
	goto tr49;
tr271:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st223;
st223:
	if ( ++p == pe )
		goto _test_eof223;
goto case; case 223:
#line 1765 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr272;
	goto tr49;
tr272:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st224;
st224:
	if ( ++p == pe )
		goto _test_eof224;
goto case; case 224:
#line 1779 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr273;
	goto tr49;
tr273:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st225;
st225:
	if ( ++p == pe )
		goto _test_eof225;
goto case; case 225:
#line 1793 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr274;
	goto tr49;
tr274:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st226;
st226:
	if ( ++p == pe )
		goto _test_eof226;
goto case; case 226:
#line 1807 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr275;
	goto tr49;
tr275:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st227;
st227:
	if ( ++p == pe )
		goto _test_eof227;
goto case; case 227:
#line 1821 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr276;
	goto tr49;
tr276:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st228;
st228:
	if ( ++p == pe )
		goto _test_eof228;
goto case; case 228:
#line 1835 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr277;
	goto tr49;
tr277:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st229;
st229:
	if ( ++p == pe )
		goto _test_eof229;
goto case; case 229:
#line 1849 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr278;
	goto tr49;
tr278:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st230;
st230:
	if ( ++p == pe )
		goto _test_eof230;
goto case; case 230:
#line 1863 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr279;
	goto tr49;
tr279:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st231;
st231:
	if ( ++p == pe )
		goto _test_eof231;
goto case; case 231:
#line 1877 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr280;
	goto tr49;
tr280:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st232;
st232:
	if ( ++p == pe )
		goto _test_eof232;
goto case; case 232:
#line 1891 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr281;
	goto tr49;
tr281:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st233;
st233:
	if ( ++p == pe )
		goto _test_eof233;
goto case; case 233:
#line 1905 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr282;
	goto tr49;
tr282:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st234;
st234:
	if ( ++p == pe )
		goto _test_eof234;
goto case; case 234:
#line 1919 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr283;
	goto tr49;
tr283:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st235;
st235:
	if ( ++p == pe )
		goto _test_eof235;
goto case; case 235:
#line 1933 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr284;
	goto tr49;
tr284:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st236;
st236:
	if ( ++p == pe )
		goto _test_eof236;
goto case; case 236:
#line 1947 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr285;
	goto tr49;
tr285:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st237;
st237:
	if ( ++p == pe )
		goto _test_eof237;
goto case; case 237:
#line 1961 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr286;
	goto tr49;
tr286:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st238;
st238:
	if ( ++p == pe )
		goto _test_eof238;
goto case; case 238:
#line 1975 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr269;
	goto tr49;
tr41:
#line 193 "sam_alignment.rl"
	{ sequence_beg = p - line.ptr; }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
goto case; case 60:
#line 1987 "sam_alignment.d"
	switch( (*p) ) {
		case 9u: goto tr101;
		case 46u: goto st60;
		case 61u: goto st60;
		default: break;
	}
	if ( (*p) > 90u ) {
		if ( 97u <= (*p) && (*p) <= 122u )
			goto st60;
	} else if ( (*p) >= 65u )
		goto st60;
	goto tr39;
tr38:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
goto case; case 61:
#line 2008 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr103;
	goto tr34;
tr103:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
goto case; case 62:
#line 2022 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr104;
	goto tr34;
tr104:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st63;
st63:
	if ( ++p == pe )
		goto _test_eof63;
goto case; case 63:
#line 2036 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr105;
	goto tr34;
tr105:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
goto case; case 64:
#line 2050 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr106;
	goto tr34;
tr106:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
goto case; case 65:
#line 2064 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr107;
	goto tr34;
tr107:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
goto case; case 66:
#line 2078 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr108;
	goto tr34;
tr108:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
goto case; case 67:
#line 2092 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr109;
	goto tr34;
tr109:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
goto case; case 68:
#line 2106 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr110;
	goto tr34;
tr110:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
goto case; case 69:
#line 2120 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr111;
	goto tr34;
tr111:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
goto case; case 70:
#line 2134 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr112;
	goto tr34;
tr112:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
goto case; case 71:
#line 2148 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr113;
	goto tr34;
tr113:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st72;
st72:
	if ( ++p == pe )
		goto _test_eof72;
goto case; case 72:
#line 2162 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr114;
	goto tr34;
tr114:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
goto case; case 73:
#line 2176 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr115;
	goto tr34;
tr115:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
goto case; case 74:
#line 2190 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr116;
	goto tr34;
tr116:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st75;
st75:
	if ( ++p == pe )
		goto _test_eof75;
goto case; case 75:
#line 2204 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr117;
	goto tr34;
tr117:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st76;
st76:
	if ( ++p == pe )
		goto _test_eof76;
goto case; case 76:
#line 2218 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr118;
	goto tr34;
tr118:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st77;
st77:
	if ( ++p == pe )
		goto _test_eof77;
goto case; case 77:
#line 2232 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr37;
	goto tr34;
tr33:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st78;
st78:
	if ( ++p == pe )
		goto _test_eof78;
goto case; case 78:
#line 2244 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr119;
	goto tr30;
tr119:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st79;
st79:
	if ( ++p == pe )
		goto _test_eof79;
goto case; case 79:
#line 2258 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr120;
	goto tr30;
tr120:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st80;
st80:
	if ( ++p == pe )
		goto _test_eof80;
goto case; case 80:
#line 2272 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr121;
	goto tr30;
tr121:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st81;
st81:
	if ( ++p == pe )
		goto _test_eof81;
goto case; case 81:
#line 2286 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr122;
	goto tr30;
tr122:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st82;
st82:
	if ( ++p == pe )
		goto _test_eof82;
goto case; case 82:
#line 2300 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr123;
	goto tr30;
tr123:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st83;
st83:
	if ( ++p == pe )
		goto _test_eof83;
goto case; case 83:
#line 2314 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr124;
	goto tr30;
tr124:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st84;
st84:
	if ( ++p == pe )
		goto _test_eof84;
goto case; case 84:
#line 2328 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr125;
	goto tr30;
tr125:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st85;
st85:
	if ( ++p == pe )
		goto _test_eof85;
goto case; case 85:
#line 2342 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr126;
	goto tr30;
tr126:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st86;
st86:
	if ( ++p == pe )
		goto _test_eof86;
goto case; case 86:
#line 2356 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr127;
	goto tr30;
tr127:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st87;
st87:
	if ( ++p == pe )
		goto _test_eof87;
goto case; case 87:
#line 2370 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr128;
	goto tr30;
tr128:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
goto case; case 88:
#line 2384 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr129;
	goto tr30;
tr129:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st89;
st89:
	if ( ++p == pe )
		goto _test_eof89;
goto case; case 89:
#line 2398 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr130;
	goto tr30;
tr130:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st90;
st90:
	if ( ++p == pe )
		goto _test_eof90;
goto case; case 90:
#line 2412 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr131;
	goto tr30;
tr131:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st91;
st91:
	if ( ++p == pe )
		goto _test_eof91;
goto case; case 91:
#line 2426 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr132;
	goto tr30;
tr132:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st92;
st92:
	if ( ++p == pe )
		goto _test_eof92;
goto case; case 92:
#line 2440 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr133;
	goto tr30;
tr133:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st93;
st93:
	if ( ++p == pe )
		goto _test_eof93;
goto case; case 93:
#line 2454 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr134;
	goto tr30;
tr134:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st94;
st94:
	if ( ++p == pe )
		goto _test_eof94;
goto case; case 94:
#line 2468 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr32;
	goto tr30;
st95:
	if ( ++p == pe )
		goto _test_eof95;
goto case; case 95:
	if ( (*p) == 9u )
		goto st14;
	goto tr24;
st96:
	if ( ++p == pe )
		goto _test_eof96;
goto case; case 96:
	if ( (*p) == 9u )
		goto tr136;
	goto tr24;
tr22:
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st97;
tr156:
#line 113 "sam_alignment.rl"
	{
        auto op = CigarOperation(cigar_op_len, cigar_op_chr);
        if (op.is_reference_consuming)
            end_pos += op.length;
        buffer.put!CigarOperation(op);
        {
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) + 1;
        }
    }
#line 28 "sam_alignment.rl"
	{ int_value = 0; }
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st97;
st97:
	if ( ++p == pe )
		goto _test_eof97;
goto case; case 97:
#line 2513 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr137;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr137:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st98;
st98:
	if ( ++p == pe )
		goto _test_eof98;
goto case; case 98:
#line 2539 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr139;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr139:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st99;
st99:
	if ( ++p == pe )
		goto _test_eof99;
goto case; case 99:
#line 2565 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr140;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr140:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st100;
st100:
	if ( ++p == pe )
		goto _test_eof100;
goto case; case 100:
#line 2591 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr141;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr141:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st101;
st101:
	if ( ++p == pe )
		goto _test_eof101;
goto case; case 101:
#line 2617 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr142;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr142:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st102;
st102:
	if ( ++p == pe )
		goto _test_eof102;
goto case; case 102:
#line 2643 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr143;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr143:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st103;
st103:
	if ( ++p == pe )
		goto _test_eof103;
goto case; case 103:
#line 2669 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr144;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr144:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st104;
st104:
	if ( ++p == pe )
		goto _test_eof104;
goto case; case 104:
#line 2695 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr145;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr145:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st105;
st105:
	if ( ++p == pe )
		goto _test_eof105;
goto case; case 105:
#line 2721 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr146;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr146:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st106;
st106:
	if ( ++p == pe )
		goto _test_eof106;
goto case; case 106:
#line 2747 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr147;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr147:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st107;
st107:
	if ( ++p == pe )
		goto _test_eof107;
goto case; case 107:
#line 2773 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr148;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr148:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st108;
st108:
	if ( ++p == pe )
		goto _test_eof108;
goto case; case 108:
#line 2799 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr149;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr149:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st109;
st109:
	if ( ++p == pe )
		goto _test_eof109;
goto case; case 109:
#line 2825 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr150;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr150:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st110;
st110:
	if ( ++p == pe )
		goto _test_eof110;
goto case; case 110:
#line 2851 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr151;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr151:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st111;
st111:
	if ( ++p == pe )
		goto _test_eof111;
goto case; case 111:
#line 2877 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr152;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr152:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st112;
st112:
	if ( ++p == pe )
		goto _test_eof112;
goto case; case 112:
#line 2903 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr153;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr153:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st113;
st113:
	if ( ++p == pe )
		goto _test_eof113;
goto case; case 113:
#line 2929 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) < 72u ) {
		if ( 48u <= (*p) && (*p) <= 57u )
			goto tr154;
	} else if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else
		goto tr138;
	goto tr20;
tr154:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st114;
st114:
	if ( ++p == pe )
		goto _test_eof114;
goto case; case 114:
#line 2955 "sam_alignment.d"
	switch( (*p) ) {
		case 61u: goto tr138;
		case 68u: goto tr138;
		case 80u: goto tr138;
		case 83u: goto tr138;
		case 88u: goto tr138;
		default: break;
	}
	if ( (*p) > 73u ) {
		if ( 77u <= (*p) && (*p) <= 78u )
			goto tr138;
	} else if ( (*p) >= 72u )
		goto tr138;
	goto tr20;
tr138:
#line 111 "sam_alignment.rl"
	{ cigar_op_len = to!uint(int_value); }
#line 112 "sam_alignment.rl"
	{ cigar_op_chr = (*p); }
	goto st115;
st115:
	if ( ++p == pe )
		goto _test_eof115;
goto case; case 115:
#line 2980 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr155;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr156;
	goto tr20;
tr19:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st116;
st116:
	if ( ++p == pe )
		goto _test_eof116;
goto case; case 116:
#line 2994 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr157;
	goto tr16;
tr157:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st117;
st117:
	if ( ++p == pe )
		goto _test_eof117;
goto case; case 117:
#line 3008 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr158;
	goto tr16;
tr158:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st118;
st118:
	if ( ++p == pe )
		goto _test_eof118;
goto case; case 118:
#line 3022 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr159;
	goto tr16;
tr159:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st119;
st119:
	if ( ++p == pe )
		goto _test_eof119;
goto case; case 119:
#line 3036 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr160;
	goto tr16;
tr160:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st120;
st120:
	if ( ++p == pe )
		goto _test_eof120;
goto case; case 120:
#line 3050 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr161;
	goto tr16;
tr161:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st121;
st121:
	if ( ++p == pe )
		goto _test_eof121;
goto case; case 121:
#line 3064 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr162;
	goto tr16;
tr162:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st122;
st122:
	if ( ++p == pe )
		goto _test_eof122;
goto case; case 122:
#line 3078 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr163;
	goto tr16;
tr163:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st123;
st123:
	if ( ++p == pe )
		goto _test_eof123;
goto case; case 123:
#line 3092 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr164;
	goto tr16;
tr164:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st124;
st124:
	if ( ++p == pe )
		goto _test_eof124;
goto case; case 124:
#line 3106 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr165;
	goto tr16;
tr165:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st125;
st125:
	if ( ++p == pe )
		goto _test_eof125;
goto case; case 125:
#line 3120 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr166;
	goto tr16;
tr166:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st126;
st126:
	if ( ++p == pe )
		goto _test_eof126;
goto case; case 126:
#line 3134 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr167;
	goto tr16;
tr167:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st127;
st127:
	if ( ++p == pe )
		goto _test_eof127;
goto case; case 127:
#line 3148 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr168;
	goto tr16;
tr168:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st128;
st128:
	if ( ++p == pe )
		goto _test_eof128;
goto case; case 128:
#line 3162 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr169;
	goto tr16;
tr169:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st129;
st129:
	if ( ++p == pe )
		goto _test_eof129;
goto case; case 129:
#line 3176 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr170;
	goto tr16;
tr170:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st130;
st130:
	if ( ++p == pe )
		goto _test_eof130;
goto case; case 130:
#line 3190 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr171;
	goto tr16;
tr171:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st131;
st131:
	if ( ++p == pe )
		goto _test_eof131;
goto case; case 131:
#line 3204 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr172;
	goto tr16;
tr172:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st132;
st132:
	if ( ++p == pe )
		goto _test_eof132;
goto case; case 132:
#line 3218 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr18;
	goto tr16;
tr15:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st133;
st133:
	if ( ++p == pe )
		goto _test_eof133;
goto case; case 133:
#line 3230 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr173;
	goto tr12;
tr173:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st134;
st134:
	if ( ++p == pe )
		goto _test_eof134;
goto case; case 134:
#line 3244 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr174;
	goto tr12;
tr174:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st135;
st135:
	if ( ++p == pe )
		goto _test_eof135;
goto case; case 135:
#line 3258 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr175;
	goto tr12;
tr175:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st136;
st136:
	if ( ++p == pe )
		goto _test_eof136;
goto case; case 136:
#line 3272 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr176;
	goto tr12;
tr176:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st137;
st137:
	if ( ++p == pe )
		goto _test_eof137;
goto case; case 137:
#line 3286 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr177;
	goto tr12;
tr177:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st138;
st138:
	if ( ++p == pe )
		goto _test_eof138;
goto case; case 138:
#line 3300 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr178;
	goto tr12;
tr178:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st139;
st139:
	if ( ++p == pe )
		goto _test_eof139;
goto case; case 139:
#line 3314 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr179;
	goto tr12;
tr179:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st140;
st140:
	if ( ++p == pe )
		goto _test_eof140;
goto case; case 140:
#line 3328 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr180;
	goto tr12;
tr180:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st141;
st141:
	if ( ++p == pe )
		goto _test_eof141;
goto case; case 141:
#line 3342 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr181;
	goto tr12;
tr181:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st142;
st142:
	if ( ++p == pe )
		goto _test_eof142;
goto case; case 142:
#line 3356 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr182;
	goto tr12;
tr182:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st143;
st143:
	if ( ++p == pe )
		goto _test_eof143;
goto case; case 143:
#line 3370 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr183;
	goto tr12;
tr183:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st144;
st144:
	if ( ++p == pe )
		goto _test_eof144;
goto case; case 144:
#line 3384 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr184;
	goto tr12;
tr184:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st145;
st145:
	if ( ++p == pe )
		goto _test_eof145;
goto case; case 145:
#line 3398 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr185;
	goto tr12;
tr185:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st146;
st146:
	if ( ++p == pe )
		goto _test_eof146;
goto case; case 146:
#line 3412 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr186;
	goto tr12;
tr186:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st147;
st147:
	if ( ++p == pe )
		goto _test_eof147;
goto case; case 147:
#line 3426 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr187;
	goto tr12;
tr187:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st148;
st148:
	if ( ++p == pe )
		goto _test_eof148;
goto case; case 148:
#line 3440 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr188;
	goto tr12;
tr188:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st149;
st149:
	if ( ++p == pe )
		goto _test_eof149;
goto case; case 149:
#line 3454 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr14;
	goto tr12;
st150:
	if ( ++p == pe )
		goto _test_eof150;
goto case; case 150:
	if ( (*p) == 9u )
		goto st6;
	goto tr7;
tr6:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st151;
st151:
	if ( ++p == pe )
		goto _test_eof151;
goto case; case 151:
#line 3473 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr190;
	goto tr3;
tr190:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st152;
st152:
	if ( ++p == pe )
		goto _test_eof152;
goto case; case 152:
#line 3487 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr191;
	goto tr3;
tr191:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st153;
st153:
	if ( ++p == pe )
		goto _test_eof153;
goto case; case 153:
#line 3501 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr192;
	goto tr3;
tr192:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st154;
st154:
	if ( ++p == pe )
		goto _test_eof154;
goto case; case 154:
#line 3515 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr193;
	goto tr3;
tr193:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st155;
st155:
	if ( ++p == pe )
		goto _test_eof155;
goto case; case 155:
#line 3529 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr194;
	goto tr3;
tr194:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st156;
st156:
	if ( ++p == pe )
		goto _test_eof156;
goto case; case 156:
#line 3543 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr195;
	goto tr3;
tr195:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st157;
st157:
	if ( ++p == pe )
		goto _test_eof157;
goto case; case 157:
#line 3557 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr196;
	goto tr3;
tr196:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st158;
st158:
	if ( ++p == pe )
		goto _test_eof158;
goto case; case 158:
#line 3571 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr197;
	goto tr3;
tr197:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st159;
st159:
	if ( ++p == pe )
		goto _test_eof159;
goto case; case 159:
#line 3585 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr198;
	goto tr3;
tr198:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st160;
st160:
	if ( ++p == pe )
		goto _test_eof160;
goto case; case 160:
#line 3599 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr199;
	goto tr3;
tr199:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st161;
st161:
	if ( ++p == pe )
		goto _test_eof161;
goto case; case 161:
#line 3613 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr200;
	goto tr3;
tr200:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st162;
st162:
	if ( ++p == pe )
		goto _test_eof162;
goto case; case 162:
#line 3627 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr201;
	goto tr3;
tr201:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st163;
st163:
	if ( ++p == pe )
		goto _test_eof163;
goto case; case 163:
#line 3641 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr202;
	goto tr3;
tr202:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st164;
st164:
	if ( ++p == pe )
		goto _test_eof164;
goto case; case 164:
#line 3655 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr203;
	goto tr3;
tr203:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st165;
st165:
	if ( ++p == pe )
		goto _test_eof165;
goto case; case 165:
#line 3669 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr204;
	goto tr3;
tr204:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st166;
st166:
	if ( ++p == pe )
		goto _test_eof166;
goto case; case 166:
#line 3683 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr205;
	goto tr3;
tr205:
#line 29 "sam_alignment.rl"
	{ int_value *= 10; int_value += (*p) - '0'; }
	goto st167;
st167:
	if ( ++p == pe )
		goto _test_eof167;
goto case; case 167:
#line 3697 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr5;
	goto tr3;
tr2:
#line 48 "sam_alignment.rl"
	{ read_name_beg = p - line.ptr; }
	goto st168;
st168:
	if ( ++p == pe )
		goto _test_eof168;
goto case; case 168:
#line 3709 "sam_alignment.d"
	if ( (*p) == 9u )
		goto tr206;
	if ( (*p) > 63u ) {
		if ( 65u <= (*p) && (*p) <= 126u )
			goto st168;
	} else if ( (*p) >= 33u )
		goto st168;
	goto tr0;
st169:
	if ( ++p == pe )
		goto _test_eof169;
goto case; case 169:
	if ( (*p) == 9u )
		goto tr209;
	goto st169;
tr209:
#line 51 "sam_alignment.rl"
	{ p--; {if (true) goto st181;} }
	goto st239;
st239:
	if ( ++p == pe )
		goto _test_eof239;
goto case; case 239:
#line 3733 "sam_alignment.d"
	goto st0;
st170:
	if ( ++p == pe )
		goto _test_eof170;
goto case; case 170:
	if ( (*p) == 9u )
		goto tr211;
	goto st170;
tr211:
#line 59 "sam_alignment.rl"
	{ p--; {if (true) goto st182;} }
	goto st240;
st240:
	if ( ++p == pe )
		goto _test_eof240;
goto case; case 240:
#line 3750 "sam_alignment.d"
	goto st0;
st171:
	if ( ++p == pe )
		goto _test_eof171;
goto case; case 171:
	if ( (*p) == 9u )
		goto tr213;
	goto st171;
tr213:
#line 68 "sam_alignment.rl"
	{ p--; {if (true) goto st183;} }
	goto st241;
st241:
	if ( ++p == pe )
		goto _test_eof241;
goto case; case 241:
#line 3767 "sam_alignment.d"
	goto st0;
st172:
	if ( ++p == pe )
		goto _test_eof172;
goto case; case 172:
	if ( (*p) == 9u )
		goto tr215;
	goto st172;
tr215:
#line 76 "sam_alignment.rl"
	{ p--; {if (true) goto st184;} }
	goto st242;
st242:
	if ( ++p == pe )
		goto _test_eof242;
goto case; case 242:
#line 3784 "sam_alignment.d"
	goto st0;
st173:
	if ( ++p == pe )
		goto _test_eof173;
goto case; case 173:
	if ( (*p) == 9u )
		goto tr217;
	goto st173;
tr217:
#line 82 "sam_alignment.rl"
	{ p--; {if (true) goto st185;} }
	goto st243;
st243:
	if ( ++p == pe )
		goto _test_eof243;
goto case; case 243:
#line 3801 "sam_alignment.d"
	goto st0;
st174:
	if ( ++p == pe )
		goto _test_eof174;
goto case; case 174:
	if ( (*p) == 9u )
		goto tr219;
	goto st174;
tr219:
#line 131 "sam_alignment.rl"
	{ p--; {if (true) goto st186;} }
	goto st244;
st244:
	if ( ++p == pe )
		goto _test_eof244;
goto case; case 244:
#line 3818 "sam_alignment.d"
	goto st0;
st175:
	if ( ++p == pe )
		goto _test_eof175;
goto case; case 175:
	if ( (*p) == 9u )
		goto tr221;
	goto st175;
tr221:
#line 163 "sam_alignment.rl"
	{ p--; {if (true) goto st187;} }
	goto st245;
st245:
	if ( ++p == pe )
		goto _test_eof245;
goto case; case 245:
#line 3835 "sam_alignment.d"
	goto st0;
st176:
	if ( ++p == pe )
		goto _test_eof176;
goto case; case 176:
	if ( (*p) == 9u )
		goto tr223;
	goto st176;
tr223:
#line 176 "sam_alignment.rl"
	{ p--; {if (true) goto st188;} }
	goto st246;
st246:
	if ( ++p == pe )
		goto _test_eof246;
goto case; case 246:
#line 3852 "sam_alignment.d"
	goto st0;
st177:
	if ( ++p == pe )
		goto _test_eof177;
goto case; case 177:
	if ( (*p) == 9u )
		goto tr225;
	goto st177;
tr225:
#line 188 "sam_alignment.rl"
	{ p--; {if (true) goto st189;} }
	goto st247;
st247:
	if ( ++p == pe )
		goto _test_eof247;
goto case; case 247:
#line 3869 "sam_alignment.d"
	goto st0;
st178:
	if ( ++p == pe )
		goto _test_eof178;
goto case; case 178:
	if ( (*p) == 9u )
		goto tr227;
	goto st178;
tr227:
#line 221 "sam_alignment.rl"
	{ p--; {if (true) goto st190;} }
	goto st248;
st248:
	if ( ++p == pe )
		goto _test_eof248;
goto case; case 248:
#line 3886 "sam_alignment.d"
	goto st0;
st179:
	if ( ++p == pe )
		goto _test_eof179;
goto case; case 179:
	if ( (*p) == 9u )
		goto tr229;
	goto st179;
tr229:
#line 251 "sam_alignment.rl"
	{ p--; {if (true) goto st251;} }
	goto st249;
st249:
	if ( ++p == pe )
		goto _test_eof249;
goto case; case 249:
#line 3903 "sam_alignment.d"
	goto st0;
st180:
	if ( ++p == pe )
		goto _test_eof180;
goto case; case 180:
	if ( (*p) == 9u )
		goto tr231;
	goto st180;
tr231:
#line 408 "sam_alignment.rl"
	{ p--; {if (true) goto st251;} }
	goto st250;
st250:
	if ( ++p == pe )
		goto _test_eof250;
goto case; case 250:
#line 3920 "sam_alignment.d"
	goto st0;
st181:
	if ( ++p == pe )
		goto _test_eof181;
goto case; case 181:
	if ( (*p) == 9u )
		goto st2;
	goto st0;
st182:
	if ( ++p == pe )
		goto _test_eof182;
goto case; case 182:
	if ( (*p) == 9u )
		goto st4;
	goto st0;
st183:
	if ( ++p == pe )
		goto _test_eof183;
goto case; case 183:
	if ( (*p) == 9u )
		goto st6;
	goto st0;
st184:
	if ( ++p == pe )
		goto _test_eof184;
goto case; case 184:
	if ( (*p) == 9u )
		goto st8;
	goto st0;
st185:
	if ( ++p == pe )
		goto _test_eof185;
goto case; case 185:
	if ( (*p) == 9u )
		goto tr235;
	goto st0;
st186:
	if ( ++p == pe )
		goto _test_eof186;
goto case; case 186:
	if ( (*p) == 9u )
		goto tr23;
	goto st0;
st187:
	if ( ++p == pe )
		goto _test_eof187;
goto case; case 187:
	if ( (*p) == 9u )
		goto st14;
	goto st0;
st188:
	if ( ++p == pe )
		goto _test_eof188;
goto case; case 188:
	if ( (*p) == 9u )
		goto st16;
	goto st0;
st189:
	if ( ++p == pe )
		goto _test_eof189;
goto case; case 189:
	if ( (*p) == 9u )
		goto st19;
	goto st0;
st190:
	if ( ++p == pe )
		goto _test_eof190;
goto case; case 190:
	if ( (*p) == 9u )
		goto st21;
	goto st0;
st251:
	if ( ++p == pe )
		goto _test_eof251;
goto case; case 251:
	if ( (*p) == 9u )
		goto st22;
	goto st0;
		default: break;
	}
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof191: cs = 191; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof192: cs = 192; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof193: cs = 193; goto _test_eof; 
	_test_eof194: cs = 194; goto _test_eof; 
	_test_eof195: cs = 195; goto _test_eof; 
	_test_eof196: cs = 196; goto _test_eof; 
	_test_eof197: cs = 197; goto _test_eof; 
	_test_eof198: cs = 198; goto _test_eof; 
	_test_eof199: cs = 199; goto _test_eof; 
	_test_eof200: cs = 200; goto _test_eof; 
	_test_eof201: cs = 201; goto _test_eof; 
	_test_eof202: cs = 202; goto _test_eof; 
	_test_eof203: cs = 203; goto _test_eof; 
	_test_eof204: cs = 204; goto _test_eof; 
	_test_eof205: cs = 205; goto _test_eof; 
	_test_eof206: cs = 206; goto _test_eof; 
	_test_eof207: cs = 207; goto _test_eof; 
	_test_eof208: cs = 208; goto _test_eof; 
	_test_eof209: cs = 209; goto _test_eof; 
	_test_eof210: cs = 210; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof211: cs = 211; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof212: cs = 212; goto _test_eof; 
	_test_eof213: cs = 213; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof214: cs = 214; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof215: cs = 215; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof216: cs = 216; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof217: cs = 217; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof218: cs = 218; goto _test_eof; 
	_test_eof219: cs = 219; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof220: cs = 220; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof221: cs = 221; goto _test_eof; 
	_test_eof222: cs = 222; goto _test_eof; 
	_test_eof223: cs = 223; goto _test_eof; 
	_test_eof224: cs = 224; goto _test_eof; 
	_test_eof225: cs = 225; goto _test_eof; 
	_test_eof226: cs = 226; goto _test_eof; 
	_test_eof227: cs = 227; goto _test_eof; 
	_test_eof228: cs = 228; goto _test_eof; 
	_test_eof229: cs = 229; goto _test_eof; 
	_test_eof230: cs = 230; goto _test_eof; 
	_test_eof231: cs = 231; goto _test_eof; 
	_test_eof232: cs = 232; goto _test_eof; 
	_test_eof233: cs = 233; goto _test_eof; 
	_test_eof234: cs = 234; goto _test_eof; 
	_test_eof235: cs = 235; goto _test_eof; 
	_test_eof236: cs = 236; goto _test_eof; 
	_test_eof237: cs = 237; goto _test_eof; 
	_test_eof238: cs = 238; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof129: cs = 129; goto _test_eof; 
	_test_eof130: cs = 130; goto _test_eof; 
	_test_eof131: cs = 131; goto _test_eof; 
	_test_eof132: cs = 132; goto _test_eof; 
	_test_eof133: cs = 133; goto _test_eof; 
	_test_eof134: cs = 134; goto _test_eof; 
	_test_eof135: cs = 135; goto _test_eof; 
	_test_eof136: cs = 136; goto _test_eof; 
	_test_eof137: cs = 137; goto _test_eof; 
	_test_eof138: cs = 138; goto _test_eof; 
	_test_eof139: cs = 139; goto _test_eof; 
	_test_eof140: cs = 140; goto _test_eof; 
	_test_eof141: cs = 141; goto _test_eof; 
	_test_eof142: cs = 142; goto _test_eof; 
	_test_eof143: cs = 143; goto _test_eof; 
	_test_eof144: cs = 144; goto _test_eof; 
	_test_eof145: cs = 145; goto _test_eof; 
	_test_eof146: cs = 146; goto _test_eof; 
	_test_eof147: cs = 147; goto _test_eof; 
	_test_eof148: cs = 148; goto _test_eof; 
	_test_eof149: cs = 149; goto _test_eof; 
	_test_eof150: cs = 150; goto _test_eof; 
	_test_eof151: cs = 151; goto _test_eof; 
	_test_eof152: cs = 152; goto _test_eof; 
	_test_eof153: cs = 153; goto _test_eof; 
	_test_eof154: cs = 154; goto _test_eof; 
	_test_eof155: cs = 155; goto _test_eof; 
	_test_eof156: cs = 156; goto _test_eof; 
	_test_eof157: cs = 157; goto _test_eof; 
	_test_eof158: cs = 158; goto _test_eof; 
	_test_eof159: cs = 159; goto _test_eof; 
	_test_eof160: cs = 160; goto _test_eof; 
	_test_eof161: cs = 161; goto _test_eof; 
	_test_eof162: cs = 162; goto _test_eof; 
	_test_eof163: cs = 163; goto _test_eof; 
	_test_eof164: cs = 164; goto _test_eof; 
	_test_eof165: cs = 165; goto _test_eof; 
	_test_eof166: cs = 166; goto _test_eof; 
	_test_eof167: cs = 167; goto _test_eof; 
	_test_eof168: cs = 168; goto _test_eof; 
	_test_eof169: cs = 169; goto _test_eof; 
	_test_eof239: cs = 239; goto _test_eof; 
	_test_eof170: cs = 170; goto _test_eof; 
	_test_eof240: cs = 240; goto _test_eof; 
	_test_eof171: cs = 171; goto _test_eof; 
	_test_eof241: cs = 241; goto _test_eof; 
	_test_eof172: cs = 172; goto _test_eof; 
	_test_eof242: cs = 242; goto _test_eof; 
	_test_eof173: cs = 173; goto _test_eof; 
	_test_eof243: cs = 243; goto _test_eof; 
	_test_eof174: cs = 174; goto _test_eof; 
	_test_eof244: cs = 244; goto _test_eof; 
	_test_eof175: cs = 175; goto _test_eof; 
	_test_eof245: cs = 245; goto _test_eof; 
	_test_eof176: cs = 176; goto _test_eof; 
	_test_eof246: cs = 246; goto _test_eof; 
	_test_eof177: cs = 177; goto _test_eof; 
	_test_eof247: cs = 247; goto _test_eof; 
	_test_eof178: cs = 178; goto _test_eof; 
	_test_eof248: cs = 248; goto _test_eof; 
	_test_eof179: cs = 179; goto _test_eof; 
	_test_eof249: cs = 249; goto _test_eof; 
	_test_eof180: cs = 180; goto _test_eof; 
	_test_eof250: cs = 250; goto _test_eof; 
	_test_eof181: cs = 181; goto _test_eof; 
	_test_eof182: cs = 182; goto _test_eof; 
	_test_eof183: cs = 183; goto _test_eof; 
	_test_eof184: cs = 184; goto _test_eof; 
	_test_eof185: cs = 185; goto _test_eof; 
	_test_eof186: cs = 186; goto _test_eof; 
	_test_eof187: cs = 187; goto _test_eof; 
	_test_eof188: cs = 188; goto _test_eof; 
	_test_eof189: cs = 189; goto _test_eof; 
	_test_eof190: cs = 190; goto _test_eof; 
	_test_eof251: cs = 251; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: 
	case 168: 
#line 50 "sam_alignment.rl"
	{ p--; {if (true) goto st169;} }
	break;
	case 2: 
	case 3: 
	case 151: 
	case 152: 
	case 153: 
	case 154: 
	case 155: 
	case 156: 
	case 157: 
	case 158: 
	case 159: 
	case 160: 
	case 161: 
	case 162: 
	case 163: 
	case 164: 
	case 165: 
	case 166: 
	case 167: 
#line 58 "sam_alignment.rl"
	{ p--; {if (true) goto st170;} }
	break;
	case 4: 
	case 5: 
	case 150: 
#line 67 "sam_alignment.rl"
	{ p--; {if (true) goto st171;} }
	break;
	case 6: 
	case 7: 
	case 133: 
	case 134: 
	case 135: 
	case 136: 
	case 137: 
	case 138: 
	case 139: 
	case 140: 
	case 141: 
	case 142: 
	case 143: 
	case 144: 
	case 145: 
	case 146: 
	case 147: 
	case 148: 
	case 149: 
#line 75 "sam_alignment.rl"
	{ p--; {if (true) goto st172;} }
	break;
	case 8: 
	case 9: 
	case 116: 
	case 117: 
	case 118: 
	case 119: 
	case 120: 
	case 121: 
	case 122: 
	case 123: 
	case 124: 
	case 125: 
	case 126: 
	case 127: 
	case 128: 
	case 129: 
	case 130: 
	case 131: 
	case 132: 
#line 81 "sam_alignment.rl"
	{ p--; {if (true) goto st173;} }
	break;
	case 10: 
	case 11: 
	case 97: 
	case 98: 
	case 99: 
	case 100: 
	case 101: 
	case 102: 
	case 103: 
	case 104: 
	case 105: 
	case 106: 
	case 107: 
	case 108: 
	case 109: 
	case 110: 
	case 111: 
	case 112: 
	case 113: 
	case 114: 
	case 115: 
#line 124 "sam_alignment.rl"
	{
        auto ptr = cast(uint*)(buffer.data.ptr + 3 * uint.sizeof);
        *ptr = (*ptr) & 0xFFFF0000;
        buffer.shrink(rollback_size);
        end_pos = pos + 1;
        p--; {if (true) goto st174;}
    }
	break;
	case 12: 
	case 13: 
	case 95: 
	case 96: 
#line 162 "sam_alignment.rl"
	{ p--; {if (true) goto st175;} }
	break;
	case 14: 
	case 15: 
	case 78: 
	case 79: 
	case 80: 
	case 81: 
	case 82: 
	case 83: 
	case 84: 
	case 85: 
	case 86: 
	case 87: 
	case 88: 
	case 89: 
	case 90: 
	case 91: 
	case 92: 
	case 93: 
	case 94: 
#line 175 "sam_alignment.rl"
	{ p--; {if (true) goto st176;} }
	break;
	case 16: 
	case 17: 
	case 18: 
	case 61: 
	case 62: 
	case 63: 
	case 64: 
	case 65: 
	case 66: 
	case 67: 
	case 68: 
	case 69: 
	case 70: 
	case 71: 
	case 72: 
	case 73: 
	case 74: 
	case 75: 
	case 76: 
	case 77: 
#line 187 "sam_alignment.rl"
	{ p--; {if (true) goto st177;} }
	break;
	case 19: 
	case 20: 
	case 60: 
#line 217 "sam_alignment.rl"
	{
        rollback_size = buffer.length;
        p--; {if (true) goto st178;}
    }
	break;
	case 21: 
#line 243 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        for (size_t i = 0; i < l_seq; ++i)
            buffer.putUnsafe!ubyte(0xFF);
        rollback_size = buffer.length;
        p--; {if (true) goto st179;}
    }
	break;
	case 25: 
	case 26: 
	case 27: 
	case 28: 
	case 29: 
	case 30: 
	case 31: 
	case 32: 
	case 33: 
	case 34: 
	case 35: 
	case 36: 
	case 37: 
	case 38: 
	case 39: 
	case 40: 
	case 41: 
	case 42: 
	case 43: 
	case 44: 
	case 45: 
	case 46: 
	case 47: 
	case 48: 
	case 49: 
	case 50: 
	case 51: 
	case 52: 
	case 53: 
	case 54: 
	case 55: 
	case 56: 
	case 57: 
	case 58: 
	case 59: 
#line 403 "sam_alignment.rl"
	{
        buffer.shrink(rollback_size);
        p--; {if (true) goto st180;}
    }
	break;
	case 192: 
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 191: 
#line 236 "sam_alignment.rl"
	{
        // '*' may correspond either to a one-base long sequence
        // or to absence of information
        if (quals_length == 1 && quals_last_char == '*' && l_seq == 0)
            buffer.shrink(rollback_size);
    }
#line 253 "sam_alignment.rl"
	{
        if (buffer.length - rollback_size != l_seq) {
            buffer.shrink(rollback_size);
            for (size_t i = 0; i < l_seq; ++i)
                buffer.putUnsafe!ubyte(0xFF);
        }
        rollback_size = buffer.length;
    }
	break;
	case 216: 
#line 326 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('Z');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 215: 
#line 337 "sam_alignment.rl"
	{
        {
        auto data = cast(ubyte[])(line[tagvalue_beg .. p - line.ptr]);
        buffer.capacity = buffer.length + 4 + data.length;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('H');
        buffer.putUnsafe(data);
        buffer.putUnsafe!ubyte(0);
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 221: 
	case 222: 
	case 223: 
	case 224: 
	case 225: 
	case 226: 
	case 227: 
	case 228: 
	case 229: 
	case 230: 
	case 231: 
	case 232: 
	case 233: 
	case 234: 
	case 235: 
	case 236: 
	case 237: 
	case 238: 
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 285 "sam_alignment.rl"
	{
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
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 193: 
	case 194: 
	case 195: 
	case 196: 
	case 197: 
	case 198: 
	case 199: 
	case 200: 
	case 201: 
	case 202: 
	case 203: 
	case 204: 
	case 205: 
	case 206: 
	case 207: 
	case 208: 
	case 209: 
	case 210: 
#line 30 "sam_alignment.rl"
	{ int_value *= current_sign; current_sign = 1; }
#line 362 "sam_alignment.rl"
	{
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
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 217: 
	case 218: 
	case 219: 
	case 220: 
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 319 "sam_alignment.rl"
	{
        buffer.capacity = buffer.length + 7;
        buffer.putUnsafe(tag_key);
        buffer.putUnsafe!char('f');
        buffer.putUnsafe!float(float_value);
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
	case 211: 
	case 212: 
	case 213: 
	case 214: 
#line 38 "sam_alignment.rl"
	{
        float_value = to!float(line[float_beg .. p - line.ptr]);
    }
#line 379 "sam_alignment.rl"
	{
        buffer.put!float(float_value);
        {
            auto ptr = cast(uint*)(buffer.data.ptr + tag_array_length_offset);
            ++*ptr;
        }
    }
#line 410 "sam_alignment.rl"
	{ rollback_size = buffer.length; }
	break;
#line 4659 "sam_alignment.d"
		default: break;
	}
	}

	_out: {}
	}

#line 481 "sam_alignment.rl"

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
