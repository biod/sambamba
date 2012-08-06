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

#line 1 "../region.rl"
module region;


#line 7 "../region.d"
static const int region_parser_start = 1;
static const int region_parser_first_final = 3;
static const int region_parser_error = 0;

static const int region_parser_en_region = 1;


#line 21 "../region.rl"


import std.conv;

struct Region {
    string reference = null;
    int beg;
    int end;
}

Region parseRegion(string str) {
    char* p = cast(char*)str.ptr;
    char* pe = p + str.length;
    char* eof = pe;
    int cs;
    long int_value;

    Region region;
    region.beg = 0;
    region.end = int.max;

    
#line 38 "../region.d"
	{
	cs = region_parser_start;
	}

#line 43 "../region.rl"
    
#line 45 "../region.d"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 1:
	if ( (*p) < 43u ) {
		if ( 33u <= (*p) && (*p) <= 41u )
			goto st3;
	} else if ( (*p) > 60u ) {
		if ( 62u <= (*p) && (*p) <= 126u )
			goto st3;
	} else
		goto st3;
	goto st0;
st0:
cs = 0;
	goto _out;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 58u )
		goto tr3;
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st3;
	goto st0;
tr3:
#line 10 "../region.rl"
	{ region.reference = str[0 .. p - str.ptr]; }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 81 "../region.d"
	if ( (*p) == 44u )
		goto tr5;
	if ( (*p) < 48u ) {
		if ( 33u <= (*p) && (*p) <= 47u )
			goto st5;
	} else if ( (*p) > 57u ) {
		if ( 58u <= (*p) && (*p) <= 126u )
			goto st5;
	} else
		goto tr5;
	goto st0;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	if ( 33u <= (*p) && (*p) <= 126u )
		goto st5;
	goto st0;
tr5:
#line 6 "../region.rl"
	{ int_value = 0; }
#line 7 "../region.rl"
	{ if ((*p) != ',') int_value *= 10, int_value += (*p) - '0'; }
	goto st6;
tr6:
#line 7 "../region.rl"
	{ if ((*p) != ',') int_value *= 10, int_value += (*p) - '0'; }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 114 "../region.d"
	switch( (*p) ) {
		case 44u: goto tr6;
		case 45u: goto tr7;
		default: break;
	}
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr6;
	goto st0;
tr7:
#line 11 "../region.rl"
	{ region.beg = to!int(int_value - 1); }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 131 "../region.d"
	if ( (*p) == 44u )
		goto tr2;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr2;
	goto st0;
tr2:
#line 6 "../region.rl"
	{ int_value = 0; }
#line 7 "../region.rl"
	{ if ((*p) != ',') int_value *= 10, int_value += (*p) - '0'; }
	goto st7;
tr8:
#line 7 "../region.rl"
	{ if ((*p) != ',') int_value *= 10, int_value += (*p) - '0'; }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 151 "../region.d"
	if ( (*p) == 44u )
		goto tr8;
	if ( 48u <= (*p) && (*p) <= 57u )
		goto tr8;
	goto st0;
		default: break;
	}
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 3: 
	case 4: 
	case 5: 
#line 10 "../region.rl"
	{ region.reference = str[0 .. p - str.ptr]; }
	break;
	case 6: 
#line 11 "../region.rl"
	{ region.beg = to!int(int_value - 1); }
	break;
	case 7: 
#line 12 "../region.rl"
	{ region.end = to!int(int_value); }
	break;
#line 184 "../region.d"
		default: break;
	}
	}

	_out: {}
	}

#line 44 "../region.rl"

    return region;
}

unittest {
    auto region1 = parseRegion("chr1:1,000-2000");
    assert(region1.reference == "chr1");
    assert(region1.beg == 999);
    assert(region1.end == 2000);

    auto region2 = parseRegion("chr2");
    assert(region2.reference == "chr2");

    auto region3 = parseRegion("chr3:1,000,000");
    assert(region3.reference == "chr3");
    assert(region3.beg == 999_999);
}
