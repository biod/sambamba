/*
    New style BGZF writer. This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/

// Based on the original by Artem Tarasov.

module bio.std.experimental.hts.bgzf_writer;

// import core.stdc.stdlib : malloc, free;
import core.memory: pureMalloc, pureFree;
import core.stdc.stdio: fopen, fread, fclose;
import std.bitmanip;
import std.conv;
import std.exception;
import std.typecons;
import std.parallelism;
import std.array;
import std.algorithm : max;
import std.stdio;
import std.typecons;

// depends on old 
import bio.core.bgzf.compress;
import bio.core.utils.roundbuf;

// import undead.stream;

import bio.std.hts.bam.constants: BGZF_MAX_BLOCK_SIZE, BGZF_BLOCK_SIZE, BGZF_EOF;
import bio.std.experimental.hts.bgzf;
import bio.std.experimental.hts.constants;

alias void delegate(ubyte[], ubyte[]) BlockWriteHandler;

/// Convenience function for Taskpool handler
Tuple!(ubyte[], ubyte[], BlockWriteHandler) bgzfCompressFunc(ubyte[] input,
                                                             int level,
                                                             ubyte[] output_buffer,
                                                             BlockWriteHandler handler)
{
  auto output = bgzfCompress(input, level, output_buffer);
  return tuple(input, output, handler);
}

/// BGZF compression - this is a port of the original that used the
/// undead.stream library.
struct BgzfWriter {

private:
  File f;
  TaskPool task_pool;

  ubyte[] buffer; // a slice into compression_buf (uncompressed data)
  ubyte[] tmp;    // a slice into compression_buf (compressed data)
  size_t current_size;
  int compression_level;

  alias Task!(bgzfCompressFunc,
              ubyte[], int, ubyte[], BlockWriteHandler) CompressionTask;
  RoundBuf!(CompressionTask*) _compression_tasks;
  ubyte[] compression_buf;

public:

  /// Create new BGZF output stream with a multi-threaded writer
  this(string fn, int _compression_level=-1) {
    f = File(fn,"wb");
    enforce1(-1 <= compression_level && compression_level <= 9,
            "BGZF compression level must be a number in interval [-1, 9]");
    size_t max_block_size = BGZF_MAX_BLOCK_SIZE;
    size_t block_size     = BGZF_BLOCK_SIZE;
    task_pool             = taskPool(),
    compression_level     = _compression_level;

    // create a ring buffer that is large enough
    size_t n_tasks = max(task_pool.size, 1) * 16;
    _compression_tasks = RoundBuf!(CompressionTask*)(n_tasks);

    // create extra block to which we can write while n_tasks are
    // executed
    auto comp_buf_size = (2 * n_tasks + 2) * max_block_size;
    auto p = cast(ubyte*)pureMalloc(comp_buf_size);
    compression_buf = p[0 .. comp_buf_size];
    buffer          = compression_buf[0 .. block_size];
    tmp             = compression_buf[max_block_size .. max_block_size * 2];
  }

  ~this() {
    close();
  }

  @disable this(this); // BgzfWriter does not have copy semantics;

  void throwBgzfException(string msg, string file = __FILE__, size_t line = __LINE__) {
    throw new BgzfException("Error writing BGZF block starting in "~f.name ~
                            " (" ~ file ~ ":" ~ to!string(line) ~ "): " ~ msg);
  }

  void enforce1(bool check, lazy string msg, string file = __FILE__, int line = __LINE__) {
    if (!check)
      throwBgzfException(msg,file,line);
  }

  void write(const void* buf, size_t size) {
    // stderr.writeln("HEY1 writing bytes ",size);
    if (size + current_size >= buffer.length) {
      size_t room;
      ubyte[] data = (cast(ubyte*)buf)[0 .. size];

      while (data.length + current_size >= buffer.length) {
        room = buffer.length - current_size;
        buffer[$ - room .. $] = data[0 .. room];
        data = data[room .. $];

        current_size = buffer.length;

        flush_block();
      }

      buffer[0 .. data.length] = data[];
      current_size = data.length;
    } else {
      buffer[current_size .. current_size + size] = (cast(ubyte*)buf)[0 .. size];
      current_size += size;
    }
    // return size;
  }

  void write(ubyte[] buf) {
    write(buf.ptr, buf.length);
  }

  void write(string s) {
    write(cast(ubyte[])s);
  }

  void write(T)(T value) { // int values
    // ubyte[T.sizeof] buf;
    ubyte[] buf = [0,0,0,0,0,0,0,0,0,0];
    assert(T.sizeof < buf.length);
    buf.write!(T,Endian.littleEndian)(value,0);
    // writeln("HEY T.sizeof: ",T.sizeof," value ",value," ",buf[0..T.sizeof]);
    write(buf[0..T.sizeof]);
  }

  /// Force flushing current block, even if it is not yet filled.
  /// Should also be used when it's not desired to have records
  /// crossing block borders.
  void flush_block() {
    if (current_size == 0)
      return;

    Tuple!(ubyte[], ubyte[], BlockWriteHandler) front_result;
    if (_compression_tasks.full) {
      front_result = _compression_tasks.front.yieldForce();
      _compression_tasks.popFront();
    }

    auto compression_task = task!bgzfCompressFunc(buffer[0 .. current_size],
                                                  compression_level, tmp,
                                                  _before_write);
    _compression_tasks.put(compression_task);
    task_pool.put(compression_task);

    size_t offset = buffer.ptr - compression_buf.ptr;
    immutable N = tmp.length;
    offset += 2 * N;
    if (offset == compression_buf.length)
      offset = 0;
    buffer = compression_buf[offset .. offset + buffer.length];
    tmp = compression_buf[offset + N .. offset + 2 * N];
    current_size = 0;

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
    // _stream.writeExact(compressed.ptr, compressed.length);
    f.rawWrite(compressed);
  }

  /// Flush all remaining BGZF blocks and underlying stream.
  void flush() {
    flush_block();

    while (!_compression_tasks.empty) {
      auto task = _compression_tasks.front;
      auto block = task.yieldForce();
      writeResult(block);
      _compression_tasks.popFront();
    }

    f.flush();
    current_size = 0;
  }

  /// Flush all remaining BGZF blocks and close source stream.
  /// Automatically adds empty block at the end, serving as indicator
  /// of end of stream. This function is automatically called on
  /// destruction.
  void close() {
    flush();
    f.rawWrite(BGZF_EOF);
    f.close();
    pureFree(compression_buf.ptr);
  }
}
