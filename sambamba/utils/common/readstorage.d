/*
    This file is part of Sambamba.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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

import core.stdc.stdlib;
import core.stdc.string;
import std.parallelism;
import std.algorithm;
import bio.core.utils.roundbuf;
import bio.core.utils.outbuffer;

import bio.bam.read;

struct ReadStorage {
    private size_t max_sz;

    this(size_t max_total_size) {
        max_sz = max_total_size;
        read_storage = cast(ubyte*)core.stdc.stdlib.malloc(max_sz);
        _reads_capa = 1024;
        auto sz = BamRead.sizeof * _reads_capa;
        _reads = cast(BamRead*)core.stdc.stdlib.malloc(sz);
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

            core.stdc.string.memcpy(read_storage + _used, read.raw_data.ptr, len);
            if (_n_reads == _reads_capa) {
                _reads_capa *= 2;
                _reads = cast(BamRead*)core.stdc.stdlib.realloc(_reads, _reads_capa * BamRead.sizeof);
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
            read_storage = cast(ubyte*)core.stdc.stdlib.realloc(read_storage, len);
            _used = len;
            read_storage[0 .. len] = read.raw_data[];
            _reads[0].raw_data = read_storage[0 .. _used];
            _reads[0].associateWithReader(read.reader);
            reads.popFront();
        }
    }

    void free() {
        core.stdc.stdlib.free(read_storage);
        core.stdc.stdlib.free(_reads);
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

class TaskWithData(alias converter, T...) {
    Task!(converter, BamRead[], OutBuffer, T)* conversion_task;
    ReadStorage input_buffer;
    OutBuffer output_buffer;

    this(size_t n=131072) {
        input_buffer = ReadStorage(n);
        output_buffer = new OutBuffer(n * 5);
    }

    final void run(TaskPool pool, T params) {
        conversion_task = task!converter(input_buffer.reads, output_buffer,
					 params);
        pool.put(conversion_task);
    }
}

void runConversion(alias converter, alias writer, R)(R reads, TaskPool pool) {
    auto n_tasks = max(pool.size, 2) * 4;
    auto tasks = RoundBuf!(TaskWithData!converter)(n_tasks);
    foreach (i; 0 .. n_tasks) {
        if (reads.empty)
            break;
        auto t = new TaskWithData!converter();
        t.input_buffer.fill(&reads);
        t.run(pool);
        tasks.put(t);
    }

    while (!tasks.empty) {
        auto t = tasks.front;
        writer(t.conversion_task.yieldForce());
        tasks.popFront();
        if (!reads.empty) {
            t.input_buffer.clear();
            t.input_buffer.fill(&reads);
            t.output_buffer.clear();
            t.run(pool);
            tasks.put(t);
        }
    }
}
