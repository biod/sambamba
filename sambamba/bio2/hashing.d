/*
    This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/

module bio2.hashing;

import std.conv;
import std.stdio;

pragma(inline):
pure nothrow @nogc ushort get16bits(char *p)
{
  ushort p0 = p[0];
  ushort p1 = p[1];
  // return p[1]*256+p[0];
  return cast(ushort)(p1*256+p0);
}

/// Paul Hsieh's fast hash (LGPL license), see http://www.azillionmonkeys.com/qed/hash.html
pure nothrow @nogc uint SuperFastHash(string str, uint hashinc = 0) {
  auto data = cast(char*)str.ptr;
  auto len  = cast(uint)str.length;

  uint hash = (hashinc > 0 ? hashinc : 0);

  if (len == 0) return 0;

  int rem = len & 3;
  len >>= 2;

  /* Main loop */
  for (;len > 0; len--) {
    hash  += get16bits(data);
    uint tmp    = (get16bits(data+2) << 11) ^ hash;
    hash   = (hash << 16) ^ tmp;
    data  += 2*ushort.sizeof;
    hash  += hash >> 11;
  }

  /* Handle end cases */
  switch (rem) {
  case 3: hash += get16bits(data);
    hash ^= hash << 16;
    hash ^= (cast(ushort)data[ushort.sizeof]) << 18;
    hash += hash >> 11;
    break;
  case 2: hash += get16bits(data);
    hash ^= hash << 11;
    hash += hash >> 17;
    break;
  case 1: hash += cast(ushort)*data;
    hash ^= hash << 10;
    hash += hash >> 1;
    break;
  default:
    break;
  }

  /* Force "avalanching" of final 127 bits */
  hash ^= hash << 3;
  hash += hash >> 5;
  hash ^= hash << 4;
  hash += hash >> 17;
  hash ^= hash << 25;
  hash += hash >> 6;

  return hash;
}


unittest {
  // Make sure it behaves the same as the original
  auto test = get16bits(cast(char *)("xy".ptr));
  assert(test == 31096, to!string(test));
  auto test2 = get16bits(cast(char *)("xyz".ptr));
  assert(test2 == 31096, to!string(test2));
  assert(get16bits(cast(char *)"xyz".ptr) == 31096);
  assert(SuperFastHash("*")==1029965590,to!string(SuperFastHash("*")));
  assert(SuperFastHash("hst")==1867544282,to!string(SuperFastHash("hst")));
  assert(SuperFastHash("hstiaashccaht")==2173265029,to!string(SuperFastHash("hstiaashccaht")));
  assert(SuperFastHash("Pjotr")==2808102289,to!string(SuperFastHash("Pjotr")));
}
