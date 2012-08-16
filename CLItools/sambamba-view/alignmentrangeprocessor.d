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
module alignmentrangeprocessor;

import bamfile;
import sam.serialize;
import bamoutput;
import jsonserialization;
import utils.format;

import std.c.stdio;

import std.exception;
import std.string;
import std.range;
import std.traits;
import std.stream;
import std.conv;
import std.range;

class ReadCounter {
    size_t number_of_reads;
    
    void process(R, SB)(R reads, SB bam) {
        number_of_reads = walkLength(reads);
    }
}

extern(C)
{
    version(Posix) {
        int isatty(int fd);

        extern(D) {
            bool isTty(shared FILE* file) {
                return isatty(fileno(cast(FILE*)file)) != 0;
            }
        }
    }

    version(Windows) {
        int _isatty(int fd);

        extern(D) {
            bool isTty(shared FILE* file) {
                return _isatty(_fileno(cast(FILE*)file)) != 0;
            }
        }
    }
}

class TextSerializer {
    this(string output_filename, bool append=false) {
        if (!isTty(stdout)) {
            // setup a buffer for stdout for faster output
            if (output_buf is null)
                output_buf = new char[1_048_576];

            setvbuf(stdout, output_buf.ptr, _IOFBF, output_buf.length);        
        }

        if (output_filename !is null)
            freopen(toStringz(output_filename), append ? "a+" : "w+", stdout);
    }

    private static __gshared char[] output_buf;

    ~this() {
        fflush(stdout);
    }
}

final class SamSerializer : TextSerializer {
    this(string filename, bool append=false) {
        super(filename, append);
    }

    void process(R, SB)(R reads, SB bam) {
        foreach (read; reads) {
            serialize(read, bam.reference_sequences, stdout);
            putcharacter(stdout, '\n');
        }
    }
}

final class BamSerializer {

    private string _output_fn;
    private int _level;

    this(string output_filename, int compression_level) {
        _output_fn = output_filename;
        _level = compression_level;
    }

    void process(R, SB)(R reads, SB bam) 
    {
        Stream output_stream;

        immutable BUFSIZE = 1_048_576;

        if (_output_fn is null) {
            output_stream = new BufferedFile(fileno(cast()stdout), FileMode.Out, BUFSIZE);
        } else {
            output_stream = new BufferedFile(_output_fn, FileMode.OutNew, BUFSIZE);
        }
        scope(exit) output_stream.close();

        writeBAM(output_stream, 
                 bam.header.text,
                 bam.reference_sequences,
                 reads,
                 _level);
    }
}

final class JsonSerializer : TextSerializer {
    this(string filename, bool append=false) {
        super(filename, append);
    }

    void process(R, SB)(R reads, SB bam) {
        foreach (read; reads) {
            jsonSerialize(read, bam.reference_sequences, stdout);
            putcharacter(stdout, '\n');
        }
    }
}
