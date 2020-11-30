/**
  (Almost) a copy-paste from contrib.undead.stream.d
*/
module bio.core.utils.switchendianness;

import core.bitop;

/***
* Switches the byte order of buffer.
* $(D size) must be even.
*/
void switchEndianness(const(void)* buffer, size_t size) 
in
{
  assert((size & 1) == 0);
}
body
{
    ubyte* startb = cast(ubyte*)buffer;
    uint* start = cast(uint*)buffer;
    switch (size) {
        case 0: break;
        case 2: {
            ubyte x = *startb;
            *startb = *(startb+1);
            *(startb+1) = x;
            break;
        }
        case 4: {
            *start = bswap(*start);
            break;
        }
        default: {
            uint* end = cast(uint*)(buffer + size - uint.sizeof);
            while (start < end) {
                uint x = bswap(*start);
                *start = bswap(*end);
                *end = x;
                ++start;
                --end;
            }
            startb = cast(ubyte*)start;
            ubyte* endb = cast(ubyte*)end;
            auto len = uint.sizeof - (startb - endb);
            if (len > 0) {
                switchEndianness(startb,len);
            }
        }
    }
}
