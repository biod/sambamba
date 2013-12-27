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
module sambamba.utils.view.alignmentrangeprocessor;

import bio.bam.reader;
import bio.bam.read;
import bio.bam.writer;
import bio.bam.thirdparty.msgpack;
import bio.core.utils.outbuffer;
import bio.core.utils.roundbuf;
import bio.core.utils.range;
import sambamba.utils.common.readstorage;

import std.stdio;
import std.exception;
import std.string;
import std.range;
import std.array;
import std.format;
import std.traits;
import std.stream : Stream, BufferedFile, FileMode;
import std.conv;
import std.algorithm;
import std.parallelism;

class ReadCounter {
    size_t number_of_reads;
    
    void process(R, SB)(R reads, SB bam) {
        number_of_reads = walkLength(reads);
    }

    enum is_serial = true;
}

version(Posix) {
    extern(C) int isatty(int fd);
}

version(Windows) {
    bool isatty(int handle) {
        return 1; // FIXME!!!
    }
}

private bool isTty(ref std.stdio.File file) @property {
    return isatty(file.fileno()) != 0;
}

template chunkToFormat(char format) {
    char[] chunkToFormat(R)(in R[] reads, OutBuffer outbuffer) {
        FormatSpec!char f = void;
        f.spec = format; // nothing else matters
        if (reads.length > 0) {
            outbuffer.capacity = reads.length * reads.front.size_in_bytes * 2;
            foreach (r; reads) {
                r.toString((const(char)[] s) { outbuffer.put(cast(ubyte[])s); }, f);
                outbuffer.put('\n');
            }
        }
        return cast(char[])(outbuffer.data);
    }
}

alias chunkToFormat!'s' chunkToSam;
alias chunkToFormat!'j' chunkToJson;

class TaskWithData(alias converter) {
    Task!(converter, BamRead[], OutBuffer)* conversion_task;
    ReadStorage input_buffer;
    OutBuffer output_buffer;

    this(size_t n=131072) {
        input_buffer = ReadStorage(n);
        output_buffer = new OutBuffer(n * 5);
    }

    final void run(TaskPool pool) {
        conversion_task = task!converter(input_buffer.reads, output_buffer);
        pool.put(conversion_task);
    }
}

void runTextConversion(alias converter, R)(R reads, TaskPool pool, File f) {
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

    auto w = f.lockingTextWriter;
    while (!tasks.empty) {
        auto t = tasks.front;
        w.put(t.conversion_task.yieldForce());
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

class TextSerializer {
    this(File f, TaskPool pool) {
        _f = f;
        _pool = pool;

        if (!_f.isTty)
            _f.setvbuf(1_024_576);
    }

    private std.stdio.File _f;
    private TaskPool _pool;
}

final class SamSerializer : TextSerializer {
    this(File f, TaskPool pool) { super(f, pool); }

    void process(R, SB)(R reads, SB bam) {
        runTextConversion!chunkToSam(reads, _pool, _f);
    }

    enum is_serial = true;
}

final class BamSerializer {

    private File _f;
    private int _level;
    private TaskPool _task_pool;
    private enum BUFSIZE = 4096;//1_048_576;


    this(File f, int compression_level, TaskPool pool) {
        _f = f;
        if (!_f.isTty)
            _f.setvbuf(BUFSIZE);
        _level = compression_level;
        _task_pool = pool;
    }

    enum is_serial = true;

    void process(R, SB)(R reads, SB bam) 
    {
        version (Posix) {
            auto handle = _f.fileno;
        }
        version (Win32) {
            import core.stdc.stdio : _fdToHandle;
            auto handle = _fdToHandle(_f.fileno);
        }
        Stream output_stream = new BufferedFile(handle, FileMode.OutNew, 
                                                BUFSIZE);
        auto writer = new BamWriter(output_stream, _level, _task_pool);
        scope(exit) writer.finish();

        writer.writeSamHeader(bam.header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);
        foreach (read; reads)
            writer.writeRecord(read);
    }
}

final class JsonSerializer : TextSerializer {
    this(File f, TaskPool pool) { super(f, pool); }

    void process(R, SB)(R reads, SB bam) {
        runTextConversion!chunkToJson(reads, _pool, _f);
    }

    enum is_serial = true;
}

final class MsgpackSerializer : TextSerializer {
    this(File f) { super(f, null); }

    void process(R, SB)(R reads, SB bam) {
        auto packer = packer(Appender!(ubyte[])());
        foreach (read; reads) {
            packer.pack(read);
            fwrite(packer.stream.data.ptr, packer.stream.data.length, ubyte.sizeof, _f.getFP());
            packer.stream.clear();
        }
    }

    enum is_serial = true;
}
