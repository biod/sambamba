/*
    This file is part of BioD.
    Copyright (C) 2013-2016    Artem Tarasov <lomereiter@gmail.com>

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
module bio.core.bgzf.inputstream;

import bio.core.bgzf.block;
import bio.core.bgzf.virtualoffset;
import bio.core.bgzf.constants;
import bio.core.bgzf.chunk;
import bio.std.hts.bam.constants;
import bio.core.utils.roundbuf;

import contrib.undead.stream;
import std.exception;
import std.conv;
import std.parallelism;
import std.array;
import std.algorithm : min, max;

/// Exception type, thrown in case of encountering corrupt BGZF blocks
class BgzfException : Exception {
    this(string msg) { super(msg); }
}

/*

Called by

randomaccessmanager.d: fillBgzfBufferFromStream(stream, true, &block, buf.ptr);
inputstream.d: auto result = fillBgzfBufferFromStream(_stream, _seekable, block, buffer,

*/

bool fillBgzfBufferFromStream(Stream stream, bool is_seekable,
                              BgzfBlock* block, ubyte* buffer,
                              size_t *number_of_bytes_read=null)
{
    if (stream.eof())
        return false;

    ulong start_offset;
    void throwBgzfException(string msg) {
        throw new BgzfException("Error reading BGZF block starting from offset " ~
                                to!string(start_offset) ~ ": " ~ msg);
    }

    if (is_seekable)
        start_offset = stream.position;

    try {
        ubyte[4] bgzf_magic = void;

        size_t bytes_read;
        while (bytes_read < 4) {
            auto buf = bgzf_magic.ptr + bytes_read;
            auto read_ = stream.read(buf[0 .. 4 - bytes_read]);
            if (read_ == 0)
                return false;
            bytes_read += read_;
        }

        if (bgzf_magic != BGZF_MAGIC) {
            throwBgzfException("wrong BGZF magic");
        }

        ushort gzip_extra_length = void;

        if (is_seekable) {
            stream.seekCur(uint.sizeof + 2 * ubyte.sizeof);
        } else {
            uint gzip_mod_time = void;
            ubyte gzip_extra_flags = void;
            ubyte gzip_os = void;
            stream.read(gzip_mod_time);
            stream.read(gzip_extra_flags);
            stream.read(gzip_os);
        }

        stream.read(gzip_extra_length);

        ushort bsize = void; // total Block SIZE minus 1
        bool found_block_size = false;

        // read extra subfields
        size_t len = 0;
        while (len < gzip_extra_length) {
            ubyte si1 = void;    // Subfield Identifier1
            ubyte si2 = void;    // Subfield Identifier2
            ushort slen = void;  // Subfield LENgth

            stream.read(si1);
            stream.read(si2);
            stream.read(slen);

            if (si1 == BAM_SI1 && si2 == BAM_SI2) {
                // found 'BC' as subfield identifier

                if (slen != 2) {
                    throwBgzfException("wrong BC subfield length: " ~
                                       to!string(slen) ~ "; expected 2");
                }

                if (found_block_size) {
                    throwBgzfException("duplicate field with block size");
                }

                // read block size
                stream.read(bsize);
                found_block_size = true;

                // skip the rest
                if (is_seekable) {
                    stream.seekCur(slen - bsize.sizeof);
                } else {
                    stream.readString(slen - bsize.sizeof);
                }
            } else {
                // this subfield has nothing to do with block size, just skip
                if (is_seekable) {
                    stream.seekCur(slen);
                } else {
                    stream.readString(slen);
                }
            }

            auto nbytes = si1.sizeof + si2.sizeof + slen.sizeof + slen;
            if (number_of_bytes_read !is null)
                *number_of_bytes_read += nbytes;
            len += nbytes;
        }

        if (len != gzip_extra_length) {
            throwBgzfException("total length of subfields in bytes (" ~
                               to!string(len) ~
                               ") is not equal to gzip_extra_length (" ~
                               to!string(gzip_extra_length) ~ ")");
        }

        if (!found_block_size) {
            throwBgzfException("block size was not found in any subfield");
        }

        // read compressed data
        auto cdata_size = bsize - gzip_extra_length - 19;
        if (cdata_size > BGZF_MAX_BLOCK_SIZE) {
            throwBgzfException("compressed data size is more than " ~
                               to!string(BGZF_MAX_BLOCK_SIZE) ~
                               " bytes, which is not allowed by " ~
                               "current BAM specification");
        }

        block.bsize = bsize;
        block.cdata_size = cast(ushort)cdata_size;

        version(extraVerbose) {
            import std.stdio;
            // stderr.writeln("[compressed] reading ", cdata_size, " bytes starting from ", start_offset);
        }
        stream.readExact(buffer, cdata_size);
        version(extraVerbose) {
            stderr.writeln("[  compressed] [write] range: ", buffer, " - ", buffer + cdata_size);
        }
        // version(extraVerbose) {stderr.writeln("[compressed] reading block crc32 and input size...");}
        stream.read(block.crc32);
        stream.read(block.input_size);

        if (number_of_bytes_read !is null)
            *number_of_bytes_read += 12 + cdata_size + block.crc32.sizeof + block.input_size.sizeof;

        // version(extraVerbose) {stderr.writeln("[compressed] read block input size: ", block.input_size);}
        block._buffer = buffer[0 .. max(block.input_size, cdata_size)];
        block.start_offset = start_offset;
        block.dirty = false;
    } catch (ReadException e) {
        throwBgzfException("stream error: " ~ e.msg);
    }

    return true;
}

///
interface BgzfBlockSupplier {
    /// Fills $(D buffer) with compressed data and points $(D block) to it.
    /// Return value is false if there is no next block.
    ///
    /// The implementation may assume that there's enough space in the buffer.
    bool getNextBgzfBlock(BgzfBlock* block, ubyte* buffer,
                          ushort* skip_start, ushort* skip_end);

    /// Total compressed size of the supplied blocks in bytes.
    /// If unknown, should return 0.
    size_t totalCompressedSize() const;
}

///
class StreamSupplier : BgzfBlockSupplier {
    private {
        Stream _stream;
        bool _seekable;
        size_t _start_offset;
        size_t _size;
        ushort _skip_start;
    }

    ///
    this(Stream stream, ushort skip_start=0) {
        _stream = stream;
        _seekable = _stream.seekable;
        _skip_start = skip_start;
        if (_seekable)
            _size = cast(size_t)(_stream.size);
    }

    ///
    bool getNextBgzfBlock(BgzfBlock* block, ubyte* buffer,
                          ushort* skip_start, ushort* skip_end) {
        auto curr_start_offset = _start_offset;

        // updates _start_offset
        auto result = fillBgzfBufferFromStream(_stream, _seekable, block, buffer,
                                               &_start_offset);
        if (!_seekable)
            block.start_offset = curr_start_offset;

        *skip_start = _skip_start;
        _skip_start = 0;
        *skip_end = 0;
        return result;
    }

    /// Stream size if available
    size_t totalCompressedSize() const {
        return _size;
    }
}

class StreamChunksSupplier : BgzfBlockSupplier {
    private {
        Stream _stream;
        Chunk[] _chunks;

        void moveToNextChunk() {
            if (_chunks.length == 0)
                return;
            size_t i = 1;
            auto beg = _chunks[0].beg;
            for ( ; i < _chunks.length; ++i)
                if (_chunks[i].beg.coffset > _chunks[0].beg.coffset)
                    break;
            _chunks = _chunks[i - 1 .. $];
            _chunks[0].beg = beg;
            _stream.seekSet(cast(size_t)_chunks[0].beg.coffset);
            version(extraVerbose) {
                import std.stdio; stderr.writeln("started processing chunk ", beg, " - ", _chunks[0].end);
            }
        }
    }

    this(Stream stream, bio.core.bgzf.chunk.Chunk[] chunks) {
        _stream = stream;
        assert(_stream.seekable);
        _chunks = chunks;
        moveToNextChunk();
    }

    ///
    bool getNextBgzfBlock(BgzfBlock* block, ubyte* buffer,
                          ushort* skip_start, ushort* skip_end)
    {
        if (_chunks.length == 0)
            return false;

        // Usually there can't be two or more chunks overlapping a
        // single block -- in such cases they are merged during
        // indexing in most implementations.
        // If this is not the case, the algorithm should still work,
        // but it might decompress the same block several times.
        //
        // On each call of this method, one of these things happen:
        // 1) We remain in the current chunk, but read next block
        // 2) We finish dealing with the current chunk, so we move to
        //    the next one. If this was the last one, false is returned.
        //
        // moveToNextChunk moves stream pointer to chunk.beg.coffset,
        // in which case skip_start should be set to chunk.beg.uoffset

        auto result = fillBgzfBufferFromStream(_stream, true, block, buffer);
        auto offset = block.start_offset;

        if (!result)
            return false;

        if (offset == _chunks[0].beg.coffset)
            *skip_start = _chunks[0].beg.uoffset; // first block in a chunk
        else
            *skip_start = 0;

        long _skip_end; // may be equal to 65536!
        if (offset == _chunks[0].end.coffset) // last block in a chunk
            _skip_end = block.input_size - _chunks[0].end.uoffset;
        else
            _skip_end = 0;

        *skip_end = cast(ushort)_skip_end;

        if (offset >= _chunks[0].end.coffset) {
            _chunks = _chunks[1 .. $];
            moveToNextChunk();
        }

        // special case: it's not actually the last block in a chunk,
        // but rather that chunk ended on the edge of two blocks
        if (block.input_size > 0 && _skip_end == block.input_size) {
            version(extraVerbose) { import std.stdio; stderr.writeln("skip_end == input size"); }
            return getNextBgzfBlock(block, buffer, skip_start, skip_end);
        }

        return true;
    }

    /// Always zero (unknown)
    size_t totalCompressedSize() const {
        return 0;
    }
}

///
// Provided an uncompressed stream block by class
class BgzfInputStream : Stream {
    private {
        BgzfBlockSupplier _supplier;
        ubyte[] _data;

        // BgzfBlockCache _cache;

        ubyte[] _read_buffer;
        VirtualOffset _current_vo;
        VirtualOffset _end_vo;

        size_t _compressed_size;

        // for estimating compression ratio
        size_t _compressed_read, _uncompressed_read;

        TaskPool _pool;
        enum _max_block_size = BGZF_MAX_BLOCK_SIZE * 2;

        alias Task!(decompressBgzfBlock, BgzfBlock) DecompressionTask;
        // DecompressionTask[] _task_buf;

        // static here means that BlockAux has no access to
        // its surrounding scope https://dlang.org/spec/struct.html
        static struct BlockAux {
          BgzfBlock block;
          ushort skip_start;
          ushort skip_end;

          DecompressionTask* task;
          // alias task this; // https://dlang.org/spec/class.html#AliasThis
        }

        RoundBuf!BlockAux _tasks = void;

        size_t _offset;

        bool fillNextBlock() {
          // Sets up a decompression task and pushes it on the roundbuf _tasks
            ubyte* p = _data.ptr + _offset;
            BlockAux b = void;
            if (_supplier.getNextBgzfBlock(&b.block, p,
                                           &b.skip_start, &b.skip_end))
            {
                if (b.block.input_size == 0) // BGZF EOF block
                    return false;

                _compressed_read += b.block.end_offset - b.block.start_offset;
                _uncompressed_read += b.block.input_size;
                version(extraVerbose) {
                    import std.stdio;
                    stderr.writeln("[creating task] ", b.block.start_offset, " / ", b.skip_start, " / ", b.skip_end);
                }

                /*
                DecompressionTask tmp = void;
                tmp = scopedTask!decompressBgzfBlock(b.block);
                auto t = _task_buf.ptr + _offset / _max_block_size;
                import core.stdc.string : memcpy;
                memcpy(t, &tmp, DecompressionTask.sizeof);
                b.task = t;
                _tasks.put(b);
                _pool.put(b.task);
                */
                // tmp = scopedTask!decompressBgzfBlock(b.block);
                auto task = task!decompressBgzfBlock(b.block);
                b.task = task;
                _tasks.put(b); // _tasks is roundbuf
                _pool.put(b.task); // _pool is thread pool
                _offset += _max_block_size;
                if (_offset == _data.length)
                    _offset = 0;
                return true;
            }
            return false;
        }

        void setupReadBuffer() {
            auto b = _tasks.front;
            auto decompressed_block = b.task.yieldForce();
            auto from = b.skip_start;
            auto to = b.block.input_size - b.skip_end;
            _read_buffer = b.block._buffer.ptr[from .. to];

            if (from == to) {
                assert(from == 0);
                setEOF();
            }

            _current_vo = VirtualOffset(b.block.start_offset, from);
            version(extraVerbose) {
                import std.stdio; stderr.writeln("[setup read buffer] ", _current_vo);
            }
            if (b.skip_end > 0)
                _end_vo = VirtualOffset(b.block.start_offset, cast(ushort)to);
            else
                _end_vo = VirtualOffset(b.block.end_offset, 0);
            _tasks.popFront();
        }

        void setEOF() {
            _current_vo = _end_vo;
            readEOF = true;
        }
    }

    this(BgzfBlockSupplier supplier,
         TaskPool pool=taskPool,
         // BgzfBlockCache cache=null,
         size_t buffer_size=0)
    {
        _supplier = supplier;
        _compressed_size = _supplier.totalCompressedSize();
        _pool = pool;
        // _cache = cache;

        // The roundbuf size (n_tasks) should be at least
        // the number of threads
        size_t n_tasks = max(pool.size, 1) * 2;
        if (buffer_size > 0)
            n_tasks = max(n_tasks, buffer_size / BGZF_MAX_BLOCK_SIZE);
        // n_tasks is 13 on my machine
        _tasks = RoundBuf!BlockAux(n_tasks);
        // _task_buf = uninitializedArray!(DecompressionTask[])(n_tasks);

        _data = uninitializedArray!(ubyte[])(n_tasks * _max_block_size);

        for (size_t i = 0; i < n_tasks; ++i)
            if (!fillNextBlock())
                break;

        if (!_tasks.empty) {
            setupReadBuffer();
        }
    }

    VirtualOffset virtualTell() const {
        return _current_vo;
    }

    override ulong seek(long offset, SeekPos whence) {
        throw new SeekException("Stream is not seekable");
    }

    override size_t writeBlock(const void* buf, size_t size) {
        throw new WriteException("Stream is not writeable");
    }

    override size_t readBlock(void* buf, size_t size) {
        version(extraVerbose) {
            import std.stdio;
            // stderr.writeln("[uncompressed] reading ", size, " bytes to address ", buf);
        }
        if (_read_buffer.length == 0) {
            assert(_tasks.empty);
            setEOF();
            return 0;
        }

        auto buffer = cast(ubyte*)buf;

        auto len = min(size, _read_buffer.length);
        buffer[0 .. len] = _read_buffer[0 .. len];
        version(extraVerbose) {
            // stderr.writeln("[uncompressed] [read] range: ", _read_buffer.ptr, " - ", _read_buffer.ptr + len);
        }
        _read_buffer = _read_buffer[len .. $];
        _current_vo = VirtualOffset(cast(ulong)_current_vo + len);

        if (_read_buffer.length == 0) {
            _current_vo = _end_vo;
            if (!_tasks.empty) {
                setupReadBuffer();
                if (!readEOF)
                    fillNextBlock();
            }
            else
                setEOF();
        }

        return len;
    }

    size_t total_compressed_size() @property const {
        return _compressed_size;
    }

    float average_compression_ratio() @property const {
        if (_compressed_read == 0)
            return 0.0;
        return cast(float)_uncompressed_read / _compressed_read;
    }
}
