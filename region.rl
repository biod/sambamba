module region;

%%{
    machine region_parser;

    action init_integer { int_value = 0; }
    action consume_next_digit { if (fc != ',') int_value *= 10, int_value += fc - '0'; }
    integer = [,0-9]+ > init_integer @consume_next_digit ;

    action set_reference { region.reference = str[0 .. p - str.ptr]; }
    action set_left_end { region.beg = to!int(int_value - 1); }
    action set_right_end { region.end = to!int(int_value); }

    reference = ([!-()+-<>-~] [!-~]*) % set_reference ;
    reference_and_left_end = reference :> ':' integer % set_left_end ;
    reference_and_both_ends = reference_and_left_end '-' integer % set_right_end ;

    region := (reference @ 0) | (reference_and_left_end @ 1) | (reference_and_both_ends @ 1);

    write data;
}%%

import std.conv;

struct Region {
    string reference;
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

    auto region3 = parseRegion("chr3:1,000,000");
    assert(region3.reference == "chr3");
    assert(region3.beg == 999_999);
}
