
#line 1 "region.rl"
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
module bio.core.region;


#line 30 "region.d"
static byte[] _region_parser_actions = [
	0, 1, 1, 1, 2, 1, 3, 1, 
	4, 2, 0, 1
];

static byte[] _region_parser_key_offsets = [
	0, 0, 6, 9, 12, 19, 21, 25
];

static char[] _region_parser_trans_keys = [
	33u, 41u, 43u, 60u, 62u, 126u, 44u, 48u, 
	57u, 58u, 33u, 126u, 44u, 33u, 47u, 48u, 
	57u, 58u, 126u, 33u, 126u, 44u, 45u, 48u, 
	57u, 44u, 48u, 57u, 0
];

static byte[] _region_parser_single_lengths = [
	0, 0, 1, 1, 1, 0, 2, 1
];

static byte[] _region_parser_range_lengths = [
	0, 3, 1, 1, 3, 1, 1, 1
];

static byte[] _region_parser_index_offsets = [
	0, 0, 4, 7, 10, 15, 17, 21
];

static byte[] _region_parser_indicies = [
	0, 0, 0, 1, 2, 2, 1, 3, 
	0, 1, 5, 4, 5, 4, 1, 4, 
	1, 6, 7, 6, 1, 8, 8, 1, 
	0
];

static byte[] _region_parser_trans_targs = [
	3, 0, 7, 4, 5, 6, 6, 2, 
	7
];

static byte[] _region_parser_trans_actions = [
	0, 0, 9, 3, 0, 9, 1, 5, 
	1
];

static byte[] _region_parser_eof_actions = [
	0, 0, 0, 3, 3, 3, 5, 7
];

static int region_parser_start = 1;
static int region_parser_first_final = 3;
static int region_parser_error = 0;

static int region_parser_en_region = 1;


#line 44 "region.rl"


import std.conv;

struct Region {
    string reference;
    uint beg;
    uint end;
}

Region parseRegion(string str) {
    char* p = cast(char*)str.ptr;
    char* pe = p + str.length;
    char* eof = pe;
    int cs;
    long uint_value;

    Region region;
    region.beg = 0;
    region.end = uint.max;

    
#line 110 "region.d"
	{
	cs = region_parser_start;
	}

#line 66 "region.rl"
    
#line 117 "region.d"
	{
	int _klen;
	uint _trans;
	byte* _acts;
	uint _nacts;
	char* _keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_keys = &_region_parser_trans_keys[_region_parser_key_offsets[cs]];
	_trans = _region_parser_index_offsets[cs];

	_klen = _region_parser_single_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += cast(uint)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _region_parser_range_lengths[cs];
	if ( _klen > 0 ) {
		char* _lower = _keys;
		char* _mid;
		char* _upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += cast(uint)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _region_parser_indicies[_trans];
	cs = _region_parser_trans_targs[_trans];

	if ( _region_parser_trans_actions[_trans] == 0 )
		goto _again;

	_acts = &_region_parser_actions[_region_parser_trans_actions[_trans]];
	_nacts = cast(uint) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 29 "region.rl"
	{ uint_value = 0; }
	break;
	case 1:
#line 30 "region.rl"
	{ if ((*p) != ',') uint_value *= 10, uint_value += (*p) - '0'; }
	break;
	case 2:
#line 33 "region.rl"
	{ region.reference = str[0 .. p - str.ptr]; }
	break;
	case 3:
#line 34 "region.rl"
	{ region.beg = to!uint(uint_value - 1); }
	break;
#line 207 "region.d"
		default: break;
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	byte* __acts = &_region_parser_actions[_region_parser_eof_actions[cs]];
	uint __nacts = cast(uint) *__acts++;
	while ( __nacts-- > 0 ) {
		switch ( *__acts++ ) {
	case 2:
#line 33 "region.rl"
	{ region.reference = str[0 .. p - str.ptr]; }
	break;
	case 3:
#line 34 "region.rl"
	{ region.beg = to!uint(uint_value - 1); }
	break;
	case 4:
#line 35 "region.rl"
	{ region.end = to!uint(uint_value); }
	break;
#line 236 "region.d"
		default: break;
		}
	}
	}

	_out: {}
	}

#line 67 "region.rl"

    return region;
}

unittest {
    auto region1 = parseRegion("chr1:1,000-2000");
    assert(region1.reference == "chr1");
    assert(region1.beg == 999);
    assert(region1.end == 2000);

    auto region2 = parseRegion("chr2");
    assert(region2.reference == "chr2");
    assert(region2.beg == 0);
    assert(region2.end == uint.max);

    auto region3 = parseRegion("chr3:1,000,000");
    assert(region3.reference == "chr3");
    assert(region3.beg == 999_999);
    assert(region3.end == uint.max);
}
