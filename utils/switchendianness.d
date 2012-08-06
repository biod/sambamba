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
/**
  Copy-pasted from std/stream.d
*/
module utils.switchendianness;

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
