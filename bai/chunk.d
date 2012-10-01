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
module bai.chunk;

import virtualoffset;

struct Chunk {
    VirtualOffset beg; /// virtual file offset of the start of the chunk
    VirtualOffset end; /// virtual file offset of the end of the chunk

    /// Compare beginnings
    int opCmp(Chunk other) const nothrow {
        if (beg < other.beg) return -1;
        if (beg > other.beg) return 1;
        if (end < other.end) return -1;
        if (end > other.end) return 1;
        return 0;
    }
}
