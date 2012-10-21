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
module region;

%%{
    machine region_parser;

    action init_integer { uint_value = 0; }
    action consume_next_digit { if (fc != ',') uint_value *= 10, uint_value += fc - '0'; }
    integer = [,0-9]+ > init_integer @consume_next_digit ;

    action set_reference { region.reference = str[0 .. p - str.ptr]; }
    action set_left_end { region.beg = to!uint(uint_value - 1); }
    action set_right_end { region.end = to!uint(uint_value); }

    reference = ([!-()+-<>-~] [!-~]*) % set_reference ;
    reference_and_left_end = reference :> ':' integer % set_left_end ;
    reference_and_both_ends = reference_and_left_end '-' integer % set_right_end ;

    region := (reference @ 0) | (reference_and_left_end @ 1) | (reference_and_both_ends @ 1);

    write data;
}%%

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

    %%write init;
    %%write exec;

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
