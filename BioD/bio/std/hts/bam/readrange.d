/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.bam.readrange;

import bio.std.hts.bam.read;
import bio.std.hts.bam.abstractreader;
import bio.std.hts.bam.reader;
import bio.core.bgzf.inputstream;
import bio.core.bgzf.virtualoffset;

import contrib.undead.stream;
import std.algorithm;
import std.system;
import std.bitmanip;

/// Read + its start/end virtual offsets
struct BamReadBlock {
    VirtualOffset start_virtual_offset; ///
    VirtualOffset end_virtual_offset; ///
    BamRead read; ///
    alias read this; ///

    ///
    BamReadBlock dup() @property const {
        return BamReadBlock(start_virtual_offset, end_virtual_offset, read.dup);
    }
}

///
mixin template withOffsets() {
    /**
        Returns: virtual offsets of beginning and end of the current read
                 plus the current read itself.
     */
    BamReadBlock front() @property {
        return BamReadBlock(_start_voffset,
                            _stream.virtualTell(),
                            _current_record);
    }

    private VirtualOffset _start_voffset;

    private void beforeNextBamReadLoad() {
        _start_voffset = _stream.virtualTell();
    }
}

///
mixin template withoutOffsets() {
    /**
        Returns: current read
     */
    ref BamRead front() @property {
        return _current_record;
    }

    private void beforeNextBamReadLoad() {}
}

/// $(D front) return type is determined by $(I IteratePolicy)
struct BamReadRange(alias IteratePolicy)
{
    /// Create new range from BgzfInputStream.
    this(BgzfInputStream stream, BamReader reader=null) {
        _stream = stream;
        _reader = reader;
        _endian_stream = new EndianStream(_stream, Endian.littleEndian);
        readNext();
    }

    ///
    bool empty() @property const {
        return _empty;
    }

    mixin IteratePolicy;

    ///
    void popFront() {
        readNext();
    }

private:
    BgzfInputStream _stream;
    EndianStream _endian_stream;

    BamReader _reader;

    BamRead _current_record;
    bool _empty = false;

    ubyte[] _buffer;

    /*
      Reads next bamRead block from stream.
     */
    void readNext() {

        // In fact, on BAM files containing a special EOF BGZF block
        // this condition will be always false!
        //
        // The reason is that we don't want to unpack next block just
        // in order to see if it's an EOF one or not.
        if (_stream.eof()) {
            _empty = true;
            return;
        }

        // In order to get the right virtual offset, we need to do it here.
        version(extraVerbose) {
            // import std.stdio; stderr.writeln("record v.o. = ", _stream.virtualTell());
        }
        beforeNextBamReadLoad();

        // Here's where _empty is really set!
        ubyte[int.sizeof] tmp = void;
        auto _read = 0;
        while (_read < int.sizeof) {
            auto _actually_read = _endian_stream.readBlock(tmp.ptr + _read, int.sizeof - _read);
            if (_actually_read == 0) {
                version(development) {
                    import std.stdio;
                    stderr.writeln("[info][bamRead range] empty, read ", _read, " bytes, expected ", int.sizeof);
                }
                _empty = true;
                return;
            }
            _read += _actually_read;
        }

        int block_size = littleEndianToNative!int(tmp);

        version(extraVerbose) {
            import std.stdio;
            stderr.writeln("[uncompressed] record size: ", block_size);
        }

        ubyte[] data = void;
        if (_reader !is null && _reader._seqprocmode) {
            if (block_size > _buffer.length)
                _buffer.length = block_size;

            data = _buffer[0 .. block_size];
        } else {
            data = allocate(block_size);
        }

        _stream.readExact(data.ptr, block_size);

        _current_record = BamRead(data);
        _current_record.associateWithReader(_reader);
    }

    private {
        ubyte[] allocate(size_t size) {
            if (_alloc_buffer_used + size > _alloc_buffer.length) {
                _alloc_buffer = uninitializedArray!(ubyte[])(max(size, 65536));
                _alloc_buffer_used = 0;
            }
            auto result = _alloc_buffer[_alloc_buffer_used .. $][0 .. size];
            _alloc_buffer_used += size;
            return result;
        }
        ubyte[] _alloc_buffer;
        size_t _alloc_buffer_used;
    }
}

/// Returns: lazy range of BamRead/BamReadBlock structs constructed from a given stream.
auto bamReadRange(alias IteratePolicy=withoutOffsets)(BgzfInputStream stream, BamReader reader) {
    return BamReadRange!IteratePolicy(stream, reader);
}
