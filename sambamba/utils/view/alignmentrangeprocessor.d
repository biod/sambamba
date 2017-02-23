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
import bio.core.utils.range;
import sambamba.utils.common.readstorage;

import std.stdio;
import std.exception;
import std.string;
import std.range;
import std.array;
import std.format;
import std.traits;
import undead.stream : Stream, BufferedFile, FileMode;
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
        auto w = _f.lockingTextWriter;
        runConversion!(chunkToSam, data => w.put(data))(reads, _pool);
    }

    enum is_serial = true;
}

import cram.writer;
final class CramSerializer {
    private CramWriter _writer;
    this(string output_fn, string ref_fn, size_t n_threads) {
        if (output_fn is null)
            output_fn = "-"; // stdout
        _writer = new CramWriter(output_fn, n_threads);
        if (ref_fn !is null)
            _writer.setFastaFilename(ref_fn);
    }

    enum is_serial = true;

    void process(R, SB)(R reads, SB bam)
    {
        scope (exit) _writer.finish();
        _writer.writeSamHeader(bam.header);
        foreach (read; reads)
            _writer.writeRecord(read);
    }
}

final class BamSerializer {

    private File _f;
    private int _level;
    private TaskPool _task_pool;
    private enum BUFSIZE = 4096;//1_048_576;

    string output_filename;

    this(string output_fn, File f, int compression_level, TaskPool pool) {
        output_filename = output_fn;
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
        writer.setFilename(output_filename);
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
        auto w = _f.lockingTextWriter;
        runConversion!(chunkToJson, data => w.put(data))(reads, _pool);
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
