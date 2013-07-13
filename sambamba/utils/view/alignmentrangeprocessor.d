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
import bio.core.utils.range;

import std.stdio;
import std.exception;
import std.string;
import std.range;
import std.format;
import std.traits;
import std.stream : Stream, BufferedFile, FileMode;
import std.conv;
import std.range;
import std.parallelism;

class ReadCounter {
    size_t number_of_reads;
    
    void process(R, SB)(R reads, SB bam) {
        number_of_reads = walkLength(reads);
    }
}

version(Posix) {
    extern(C) int isatty(int fd);
}

version(Windows) {
    extern(C) int _isatty(int fd);
    alias _isatty isatty;
}

private bool isTty(ref std.stdio.File file) @property {
    return isatty(file.fileno()) != 0;
}

template chunkToFormat(string format) {
    enum spec = FormatSpec!char(format);
    char[] chunkToFormat(R)(in R[] reads) {
        auto buf = Appender!(char[])();
        if (reads.length == 0)
            return buf.data;
        buf.reserve(reads.length * reads.front.size_in_bytes * 2);
        foreach (read; reads) {
            read.toString((const(char)[] s) { buf.put(s); }, spec);
            buf.put('\n');
        }
        return buf.data;
    }
}

alias chunkToFormat!"%s" chunkToSam;
alias chunkToFormat!"%j" chunkToJson;

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
        auto read_batches = chunked(reads, 1024);
        auto sam_chunks = _pool.map!chunkToSam(read_batches, 16);
        auto w = _f.lockingTextWriter;
        foreach (chunk; sam_chunks)
            w.put(chunk);
    }
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

    void process(R, SB)(R reads, SB bam) 
    {
        Stream output_stream = new BufferedFile(_f.fileno(), FileMode.OutNew, 
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
        auto read_batches = chunked(reads, 1024);
        auto json_chunks = _pool.map!chunkToJson(read_batches, 16);
        auto w = _f.lockingTextWriter;
        foreach (chunk; json_chunks)
            w.put(chunk);
    }
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
}
