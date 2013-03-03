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
module sambamba.utils.view.alignmentrangeprocessor;

import bio.bam.reader;
import bio.bam.writer;
import bio.bam.thirdparty.msgpack;

import std.stdio;
import std.exception;
import std.string;
import std.range;
import std.format;
import std.traits;
import std.stream : Stream, BufferedFile, FileMode;
import std.conv;
import std.range;

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

class TextSerializer {
    this(File f) {
        _f = f;

        if (!_f.isTty)
            _f.setvbuf(1_024_576);
    }

    private std.stdio.File _f;
}

final class SamSerializer : TextSerializer {
    this(File f) { super(f); }
    void process(R, SB)(R reads, SB bam) {
        _f.lockingTextWriter.formattedWrite("%(%s\n%)\n", reads);
    }
}

final class BamSerializer {

    private File _f;
    private int _level;

    this(File f, int compression_level) {
        _f = f;
        _level = compression_level;
    }

    void process(R, SB)(R reads, SB bam) 
    {
        immutable BUFSIZE = 1_048_576;
        Stream output_stream = new BufferedFile(_f.fileno(), FileMode.OutNew, 
                                                BUFSIZE);
        auto writer = new BamWriter(output_stream, _level);
        scope(exit) writer.finish();

        writer.writeSamHeader(bam.header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);
        foreach (read; reads)
            writer.writeRecord(read);
    }
}

final class JsonSerializer : TextSerializer {
    this(File f) { super(f); }
    void process(R, SB)(R reads, SB bam) {
        _f.lockingTextWriter.formattedWrite("%(%j\n%)\n", reads);
    }
}

final class MsgpackSerializer : TextSerializer {
    this(File f) { super(f); }
    void process(R, SB)(R reads, SB bam) {
        auto packer = packer(Appender!(ubyte[])());
        foreach (read; reads) {
            packer.pack(read);
            fwrite(packer.stream.data.ptr, packer.stream.data.length, ubyte.sizeof, _f.getFP());
            packer.stream.clear();
        }
    }
}
