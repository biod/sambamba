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
