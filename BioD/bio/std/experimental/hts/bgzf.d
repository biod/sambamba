/*
    This file is part of Sambamba.
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

module bio.std.experimental.hts.bgzf;

/**

   This is a new version of sambamba bgzf (under development). Bgzf is
   a blocked version of the ubiquous gzip format. By making it blocked
   it allows for seeking in gz files. Note that without seeking it
   can take standard gz files too.

   The new version is a prototype for new sambamba architecture using
   canonical D language features, including immutable and improved
   laziness and a more functional programming style. It should provide
   improved performance and minimize RAM use, as well as better
   composability.

   Authors: Pjotr Prins

 */

import core.stdc.string : memcpy;

import std.bitmanip;
import std.conv;
import std.exception;
import std.file;
import std.stdio;
import std.typecons;
import std.zlib : calc_crc32 = crc32, ZlibException;

import bio.std.hts.bam.constants;
import bio.core.bgzf.block;
import bio.core.bgzf.constants;
import bio.core.utils.zlib : inflateInit2, inflate, inflateEnd, Z_OK, Z_FINISH, Z_STREAM_END;

import bio.std.experimental.hts.constants;

class BgzfException : Exception {
    this(string msg) { super(msg); }
}

alias Nullable!ulong FilePos;
alias immutable(uint) CRC32;

alias BGZF_MAX_BLOCK_SIZE BLOCK_SIZE;

alias ubyte[BLOCK_SIZE] BlockBuffer;

@property  ubyte read_ubyte(File f) {
    ubyte[1] ubyte1; // read buffer
    immutable ubyte[1] buf = f.rawRead(ubyte1);
    return buf[0];
}

@property ushort read_ushort(File f) {
    ubyte[2] ubyte2; // read buffer
    immutable ubyte[2] buf = f.rawRead(ubyte2);
    return littleEndianToNative!ushort(buf);
  }
@property uint read_uint(File f) {
    ubyte[4] ubyte4; // read buffer
    immutable ubyte[4] buf = f.rawRead(ubyte4);
    return littleEndianToNative!uint(buf);
  }

/**
   Uncompress a zlib buffer (without header)
*/
ubyte[] deflate(ubyte[] uncompressed_buf, const ubyte[] compressed_buf, size_t uncompressed_size, CRC32 crc32) {
  assert(uncompressed_buf.length == BLOCK_SIZE);
  bio.core.utils.zlib.z_stream zs;
  zs.next_in = cast(typeof(zs.next_in))compressed_buf;
  zs.avail_in = to!uint(compressed_buf.length);

  auto err = inflateInit2(&zs, /* winbits = */-15);
  if (err != Z_OK) throw new ZlibException(err);

  zs.next_out = cast(typeof(zs.next_out))uncompressed_buf.ptr;
  zs.avail_out = cast(int)uncompressed_buf.length;

  scope(exit) { inflateEnd(&zs); }
  err = inflate(&zs, Z_FINISH);
  if (err != Z_STREAM_END) throw new ZlibException(err);

  assert(zs.total_out == uncompressed_size);
  uncompressed_buf.length = uncompressed_size;
  assert(crc32 == calc_crc32(0, uncompressed_buf[]));

  return uncompressed_buf;
}

/**
    BgzfReader is designed to run on a single thread. All it does is
    fetch block headers and data, so the thread should easily keep up
    with IO. All data processing is happening lazily in other threads.
*/
struct BgzfReader {
  File f;
  FilePos report_fpos; // for error handler - assumes one thread!

  this(string fn) {
    enforce(fn.isFile);
    f = File(fn,"r");
  }

  @disable this(this); // BgzfReader does not have copy semantics;

  void throwBgzfException(string msg, string file = __FILE__, size_t line = __LINE__) {
    throw new BgzfException("Error reading BGZF block starting in "~f.name ~" @ " ~
                            to!string(report_fpos) ~ " (" ~ file ~ ":" ~ to!string(line) ~ "): " ~ msg);
  }

  void enforce1(bool check, lazy string msg, string file = __FILE__, int line = __LINE__) {
    if (!check)
      throwBgzfException(msg,file,line);
  }

  /**
      Reads the block header and returns the contained compressed data
      size with the file pointer positioned at the associated
      compressed data.
  */
  size_t read_block_header() {
    ubyte[4] ubyte4;
    auto magic = f.rawRead(ubyte4);
    enforce1(magic.length == 4, "Premature end of file");
    enforce1(magic[0..4] == BGZF_MAGIC,"Invalid file format: expected bgzf magic number");
    ubyte[uint.sizeof + 2 * ubyte.sizeof] skip;
    f.rawRead(skip); // skip gzip info
    ushort gzip_extra_length = f.read_ushort();
    immutable fpos1 = f.tell;
    size_t bsize = 0;
    while (f.tell < fpos1 + gzip_extra_length) {
      immutable subfield_id1 = f.read_ubyte();
      immutable subfield_id2 = f.read_ubyte();
      immutable subfield_len = f.read_ushort();
      if (subfield_id1 == BAM_SI1 && subfield_id2 == BAM_SI2) {
        // BC identifier
        enforce(gzip_extra_length == 6);
        // FIXME: always picks first BC block
        bsize = 1+f.read_ushort(); // BLOCK size
        enforce1(subfield_len == 2, "BC subfield len should be 2");
        break;
      }
      else {
        f.seek(subfield_len,SEEK_CUR);
      }
      enforce1(bsize!=0,"block size not found");
      f.seek(fpos1+gzip_extra_length); // skip any extra subfields - note we don't check for second BC
    }
    immutable compressed_size = bsize - 1 - gzip_extra_length - 19;
    enforce1(compressed_size <= BLOCK_SIZE, "compressed size larger than allowed");

    // stderr.writeln("[compressed] size ", compressed_size, " bytes starting block @ ", report_fpos);
    return compressed_size;
  }

  /**
     Fetch the compressed data part of the block and return it with
     the uncompressed size and CRC32. The file pointer is assumed to
     be at the start of the compressed data and will be at the end of
     that section after.
  */
  Tuple!(ubyte[],immutable(uint),CRC32) read_compressed_data(ubyte[] buffer) {
    auto compressed_buf = f.rawRead(buffer);

    immutable CRC32 crc32 = f.read_uint();
    immutable uncompressed_size = f.read_uint();
    // stderr.writeln("[uncompressed] size ",uncompressed_size);
    return tuple(compressed_buf,uncompressed_size,crc32);
  }

  /**
   * Returns new tuple of the new file position, the compressed buffer and
   * the CRC32 o the uncompressed data. file pos is NULL when done
   */
  Tuple!(FilePos,ubyte[],size_t,CRC32) read_compressed_block(FilePos fpos, ubyte[] buffer) {
    immutable start_offset = fpos;
    try {
      if (fpos.isNull) throwBgzfException("Trying to read past eof");
      report_fpos = fpos;
      f.seek(fpos);
      immutable compressed_size = read_block_header();
      auto ret = read_compressed_data(buffer[0..compressed_size]);
      auto compressed_buf = ret[0];
      immutable uncompressed_size = ret[1];
      immutable crc32 = ret[2];

      if (uncompressed_size == 0) {
        // check for eof marker, rereading block header
        auto lastpos = f.tell();
        f.seek(start_offset);
        ubyte[BGZF_EOF.length] buf;
        f.rawRead(buf);
        f.seek(lastpos);
        if (buf == BGZF_EOF)
          return tuple(FilePos(),compressed_buf,cast(size_t)0,crc32); // sets fpos to null
      }
      return tuple(FilePos(f.tell()),compressed_buf,cast(size_t)uncompressed_size,crc32);
    } catch (Exception e) { throwBgzfException(e.msg,e.file,e.line); }
    assert(0); // never reached
  }
}

/**
   Simple block iterator
*/
struct BgzfBlocks {
  BgzfReader bgzf;

  this(string fn) {
    bgzf = BgzfReader(fn);
  }

  @disable this(this); // disable copy semantics;

  int opApply(scope int delegate(ubyte[]) dg) {
    FilePos fpos = 0;

    try {
      while (!fpos.isNull) {
        BlockBuffer stack_buffer;
        auto res = bgzf.read_compressed_block(fpos,stack_buffer);
        fpos = res[0]; // point fpos to next block
        if (fpos.isNull) break;

        auto compressed_buf = res[1]; // same as stack_buffer
        auto uncompressed_size = res[2];
        auto crc32 = res[3];
        BlockBuffer uncompressed_buf;
        // call delegated function with new block
        dg(deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32));
      }
    } catch (Exception e) { bgzf.throwBgzfException(e.msg,e.file,e.line); }
    return 0;
  }
}


Tuple!(size_t,FilePos) read_blockx(ref BgzfReader bgzf, FilePos fpos, ref ubyte[] uncompressed_buf) {
  BlockBuffer compressed_buf;
  auto res = bgzf.read_compressed_block(fpos,compressed_buf);
  fpos = res[0]; // point fpos to next block
  if (fpos.isNull) return tuple(cast(size_t)0,fpos);
  auto data = res[1];

  assert(data.ptr == compressed_buf.ptr);
  size_t uncompressed_size = res[2];
  auto crc32 = res[3];
  deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32);
  return tuple(uncompressed_size,fpos);
}

import std.parallelism;

int kick_off_reading_block_ahead(ubyte[] uncompressed_buf, ubyte[] compressed_buf, size_t uncompressed_size, CRC32 crc32) {
  // writeln("HEY " ~ to!string(uncompressed_size));
  deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32);
  return -1;
}

/**
*/
struct BlockReadAhead {
  bool task_running = false, we_have_a_task = false;
  Task!(kick_off_reading_block_ahead, ubyte[], ubyte[], size_t, CRC32)* t;
  FilePos fpos2;
  size_t uncompressed_size2 = 0;
  BlockBuffer compressed_buf2;
  BlockBuffer uncompressed_buf2;

  private void read_next_block() {
  }

  private void add_deflate_task() {
  }

  private void copy_deflated_buffer() {
  }

  void setup_block_reader() {
    read_next_block();
    add_deflate_task();
    throw new Exception("NYI");
  }

  Tuple!(size_t,FilePos) read_block(ref BgzfReader bgzf, FilePos fpos, ref ubyte[] uncompressed_buf) {
    assert(we_have_a_task);
    copy_deflated_buffer();
    read_next_block();
    add_deflate_task();
    // return

    if (task_running) {
      int res = t.yieldForce;
      // writeln(res);
      task_running = false;
      memcpy(uncompressed_buf.ptr,compressed_buf2.ptr,uncompressed_size2);
      return tuple(uncompressed_size2, fpos2);
    }
    else {
      BlockBuffer compressed_buf;
      auto res = bgzf.read_compressed_block(fpos,compressed_buf);
      fpos = res[0]; // point fpos to next block
      if (fpos.isNull) return tuple(cast(size_t)0,fpos);
      auto data = res[1];
      assert(data.ptr == compressed_buf.ptr);
      size_t uncompressed_size = res[2];
      auto crc32 = res[3];

      deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32);

      // now set up a new buffer
      auto res2 = bgzf.read_compressed_block(fpos,compressed_buf2);
      fpos2 = res[0]; // point fpos to next block
      if (!fpos2.isNull) {
        auto data2 = res2[1];
        uncompressed_size2 = res2[2];
        t = task!kick_off_reading_block_ahead(cast(ubyte[])uncompressed_buf2,cast(ubyte[])compressed_buf2,uncompressed_size2,res2[3]);
        t.executeInNewThread();
        task_running = true;
      }
      return tuple(uncompressed_size,fpos);
    }
  }
}

/**
*/
struct BlockReadUnbuffered {

  void setup_block_reader() {
  }

  Tuple!(ubyte[], size_t, FilePos) read_block(ref BgzfReader bgzf, in FilePos fpos, ubyte[] uncompressed_buf) {
    BlockBuffer compressed_buf;
    auto res = bgzf.read_compressed_block(fpos,compressed_buf);
    auto fpos2 = res[0]; // point fpos to next block
    if (fpos.isNull) return tuple(uncompressed_buf,cast(size_t)0,fpos2);
    auto data = res[1];
    assert(data.ptr == compressed_buf.ptr);
    size_t uncompressed_size = res[2];
    auto crc32 = res[3];

    auto buf = deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32);
    assert(buf.ptr == uncompressed_buf.ptr);
    return tuple(uncompressed_buf,uncompressed_size,fpos2);
  }
}

/**
   Streams bgzf data and fetch items by unit or buffer. These can go beyond
   the size of individual blocks(!)
*/

struct BgzfStream {
  BgzfReader bgzf;
  FilePos fpos;             // track file position
  ubyte[] uncompressed_buf; // current data buffer
  size_t uncompressed_size; // current data buffer size
  Nullable!int block_pos;   // position in block
  BlockReadUnbuffered blockread;

  this(string fn) {
    bgzf = BgzfReader(fn);
    uncompressed_buf = new ubyte[BLOCK_SIZE];
    fpos = 0;
  }

  @disable this(this); // disable copy semantics;

  @property bool eof() {
    return fpos.isNull;
  }

  /**
     Fetch data into buffer. The size of the buffer can be larger than
     one or more multiple blocks
  */
  ubyte[] fetch(ubyte[] buffer) {
    if (block_pos.isNull) {
      blockread.setup_block_reader();
      auto res = blockread.read_block(bgzf,fpos,uncompressed_buf); // read first block
      assert(res[0].ptr == uncompressed_buf.ptr);
      uncompressed_size = res[1];
      fpos = res[2];
      block_pos = 0;
    }

    immutable buffer_length = buffer.length;
    size_t buffer_pos = 0;
    size_t remaining = buffer_length;

    while (remaining > 0) {
      if (block_pos + remaining < uncompressed_size) {
        // full copy
        assert(buffer_pos + remaining == buffer_length);
        memcpy(buffer[buffer_pos..buffer_pos+remaining].ptr,uncompressed_buf[block_pos..block_pos+remaining].ptr,remaining);
        block_pos += remaining;
        remaining = 0;
      }
      else {
        // read tail of buffer
        immutable tail = uncompressed_size - block_pos;
        memcpy(buffer[buffer_pos..buffer_pos+tail].ptr,uncompressed_buf[block_pos..uncompressed_size].ptr,tail);
        buffer_pos += tail;
        remaining -= tail;
        auto res = blockread.read_block(bgzf,fpos,uncompressed_buf);
        assert(res[0].ptr == uncompressed_buf.ptr);
        uncompressed_size = res[1];
        fpos = res[2];
        block_pos = 0;
      }
    }
    return buffer;
  }

  int read(T)() { // for integers
    ubyte[T.sizeof] buf;
    auto b = fetch(buf);
    return b.read!(T,Endian.littleEndian)();
  }

  string read(T)(size_t len) {
    ubyte[] buf = new ubyte[len]; // heap allocation
    fetch(buf);
    return cast(T)buf;
  }

  T[] read(T)(T[] buffer) { return cast(T[])fetch(cast(ubyte[])buffer); };
}
