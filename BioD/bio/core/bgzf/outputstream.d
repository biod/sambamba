/*
    This file is part of BioD.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
module bio.core.bgzf.outputstream;

import bio.core.bgzf.constants;
import bio.core.bgzf.compress;

import bio.core.utils.roundbuf;

import contrib.undead.stream;
import std.exception;
import std.parallelism;
import std.array;
import std.algorithm : max;
import std.typecons;
import core.stdc.stdlib;

alias void delegate(ubyte[], ubyte[]) BlockWriteHandler;

Tuple!(ubyte[], ubyte[], BlockWriteHandler)
bgzfCompressFunc(ubyte[] input, int level, ubyte[] output_buffer,
		 BlockWriteHandler handler)
{
    auto output = bgzfCompress(input, level, output_buffer);
    return tuple(input, output, handler);
}

/// Class for BGZF compression
class BgzfOutputStream : Stream {

    private {
        Stream _stream = void;
        TaskPool _task_pool = void;

        ubyte[] _buffer; // a slice into _compression_buffer (uncompressed data)
        ubyte[] _tmp;    // a slice into _compression_buffer (compressed data)
        size_t _current_size;

        int _compression_level;

        alias Task!(bgzfCompressFunc,
		    ubyte[], int, ubyte[], BlockWriteHandler) CompressionTask;
        RoundBuf!(CompressionTask*) _compression_tasks;
        ubyte[] _compression_buffer;
    }

    /// Create new BGZF output stream which will use
    /// provided $(D task_pool) to do multithreaded compression.
    this(Stream output_stream,
         int compression_level=-1,
         TaskPool task_pool=taskPool,
         size_t buffer_size=0,
         size_t max_block_size=BGZF_MAX_BLOCK_SIZE,
         size_t block_size=BGZF_BLOCK_SIZE)
    {
        enforce(-1 <= compression_level && compression_level <= 9,
                "Compression level must be a number in interval [-1, 9]");
        _stream = output_stream;
        _task_pool = task_pool;
        _compression_level = compression_level;

        size_t n_tasks = max(task_pool.size, 1) * 16;
        if (buffer_size > 0) {
            n_tasks = max(n_tasks, buffer_size / max_block_size);
        }
        _compression_tasks = RoundBuf!(CompressionTask*)(n_tasks);

        // 1 extra block to which we can write while n_tasks are executed
        auto comp_buf_size = (2 * n_tasks + 2) * max_block_size;
        auto p = cast(ubyte*)core.stdc.stdlib.malloc(comp_buf_size);
        _compression_buffer = p[0 .. comp_buf_size];
        _buffer = _compression_buffer[0 .. block_size];
        _tmp = _compression_buffer[max_block_size .. max_block_size * 2];

        readable = false;
        writeable = true;
        seekable = false;
    }

    override size_t readBlock(void* buffer, size_t size) {
        throw new ReadException("Stream is not readable");
    }

    override ulong seek(long offset, SeekPos whence) {
        throw new SeekException("Stream is not seekable");
    }

    override size_t writeBlock(const void* buf, size_t size) {
        if (size + _current_size >= _buffer.length) {
            size_t room;
            ubyte[] data = (cast(ubyte*)buf)[0 .. size];

            while (data.length + _current_size >= _buffer.length) {
                room = _buffer.length - _current_size;
                _buffer[$ - room .. $] = data[0 .. room];
                data = data[room .. $];

                _current_size = _buffer.length;

                flushCurrentBlock();
            }

            _buffer[0 .. data.length] = data[];
            _current_size = data.length;
        } else {
            _buffer[_current_size .. _current_size + size] = (cast(ubyte*)buf)[0 .. size];
            _current_size += size;
        }

        return size;
    }

    /// Force flushing current block, even if it is not yet filled.
    /// Should be used when it's not desired to have records crossing block borders.
    void flushCurrentBlock() {

        if (_current_size == 0)
            return;

        Tuple!(ubyte[], ubyte[], BlockWriteHandler) front_result;
        if (_compression_tasks.full) {
            front_result = _compression_tasks.front.yieldForce();
            _compression_tasks.popFront();
        }

        auto compression_task = task!bgzfCompressFunc(_buffer[0 .. _current_size],
						      _compression_level, _tmp,
						      _before_write);
        _compression_tasks.put(compression_task);
        _task_pool.put(compression_task);

        size_t offset = _buffer.ptr - _compression_buffer.ptr;
        immutable N = _tmp.length;
        offset += 2 * N;
        if (offset == _compression_buffer.length)
            offset = 0;
        _buffer = _compression_buffer[offset .. offset + _buffer.length];
        _tmp = _compression_buffer[offset + N .. offset + 2 * N];
        _current_size = 0;

        if (front_result[0] !is null)
	    writeResult(front_result);

        while (!_compression_tasks.empty) {
            auto task = _compression_tasks.front;
            if (!task.done())
                break;
            auto result = task.yieldForce();
	    writeResult(result);
            _compression_tasks.popFront();
        }
    }

    private void delegate(ubyte[], ubyte[]) _before_write;
    void setWriteHandler(void delegate(ubyte[], ubyte[]) handler) {
	_before_write = handler;
    }

    private void writeResult(Tuple!(ubyte[], ubyte[], BlockWriteHandler) block) {
	auto uncompressed = block[0];
	auto compressed = block[1];
	auto handler = block[2];
	if (handler) {// write handler enabled
	    handler(uncompressed, compressed);
	}
	_stream.writeExact(compressed.ptr, compressed.length);
    }

    /// Flush all remaining BGZF blocks and underlying stream.
    override void flush() {
        flushCurrentBlock();

	while (!_compression_tasks.empty) {
            auto task = _compression_tasks.front;
            auto block = task.yieldForce();
            writeResult(block);
            _compression_tasks.popFront();
        }

        _stream.flush();
        _current_size = 0;
    }

    /// Flush all remaining BGZF blocks and close source stream.
    /// Automatically adds empty block at the end, serving as
    /// indicator of end of stream.
    override void close() {
        flush();

        addEofBlock();

        _stream.close();

        writeable = false;
        core.stdc.stdlib.free(_compression_buffer.ptr);
    }

    /// Adds EOF block. This function is called in close() method.
    void addEofBlock() {
        _stream.writeExact(BGZF_EOF.ptr, BGZF_EOF.length);
    }
}

unittest {
    import bio.core.bgzf.inputstream;

    import std.array, std.range, std.stdio;

    char[] data = "my very l" ~ array(repeat('o', 1000000)) ~ "ng string";

    foreach (level; [-1, 0, 1]) {
        auto output_stream = new MemoryStream();
        auto bgzf_output_stream = new BgzfOutputStream(output_stream, 1);
        bgzf_output_stream.writeExact(data.ptr, data.length);
        bgzf_output_stream.close();

        auto input_stream = new MemoryStream(output_stream.data);
        input_stream.seekSet(0);

        auto block_supplier = new StreamSupplier(input_stream);
        auto bgzf_input_stream = new BgzfInputStream(block_supplier);
        char[] uncompressed_data = new char[data.length];
        bgzf_input_stream.readExact(uncompressed_data.ptr, data.length);
        assert(uncompressed_data == data);
    }
}
