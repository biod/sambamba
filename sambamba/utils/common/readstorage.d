/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.utils.common.readstorage;

import std.c.stdlib;
import std.c.string;

import bio.bam.read;

struct ReadStorage {
    private size_t max_sz;

    this(size_t max_total_size) {
        max_sz = max_total_size;
        read_storage = cast(ubyte*)std.c.stdlib.malloc(max_sz);
        _reads_capa = 1024;
        auto sz = BamRead.sizeof * _reads_capa;
        _reads = cast(BamRead*)std.c.stdlib.malloc(sz);
    }

    void clear() {
        _used = 0;
        _n_reads = 0;
    }

    void fill(R)(R* reads) {
        while (!reads.empty) {
            auto read = reads.front;
            auto len = read.raw_data.length;
            if (len + _used > max_sz)
                break;

            std.c.string.memcpy(read_storage + _used, read.raw_data.ptr, len);
            if (_n_reads == _reads_capa) {
                _reads_capa *= 2;
                _reads = cast(BamRead*)std.c.stdlib.realloc(_reads, _reads_capa * BamRead.sizeof);
            }
            _reads[_n_reads].raw_data = read_storage[_used .. _used + len];
            _reads[_n_reads].associateWithReader(read.reader);

            _n_reads += 1;
            _used += len;

            reads.popFront();
        }

        if (_n_reads == 0) {
            auto read = reads.front;
            auto len = read.raw_data.length;
            assert(len > max_sz);
            _n_reads = 1;
            read_storage = cast(ubyte*)std.c.stdlib.realloc(read_storage, len);
            _used = len;
            read_storage[0 .. len] = read.raw_data[];
            _reads[0].raw_data = read_storage[0 .. _used];
            _reads[0].associateWithReader(read.reader);
            reads.popFront();
        }
    }

    void free() {
        std.c.stdlib.free(read_storage);
        std.c.stdlib.free(_reads);
    }

    BamRead[] reads() @property { return _reads[0 .. _n_reads]; }

    private {
        BamRead* _reads;
        size_t _reads_capa;
        size_t _n_reads;

        ubyte* read_storage;
        size_t _used;
    }
}
