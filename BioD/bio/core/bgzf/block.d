/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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
module bio.core.bgzf.block;

import bio.std.hts.bam.constants;
// import bio.core.utils.memoize;
import bio.core.utils.zlib;

import std.array;
import std.conv;
import std.algorithm;
import std.exception;

/**
  Structure representing BGZF block.
  In general, users shouldn't use it, as it is EXTREMELY low-level.

  Note it is a struct that has support for comparison based
  on its crc32 value.
 */
struct BgzfBlock {
    // field types are as in the SAM/BAM specification
    // ushort ~ uint16_t, char ~ uint8_t, uint ~ uint32_t

    public ulong start_offset; /// start offset in the file, in bytes

    /// end offset in the file, in bytes
    public ulong end_offset() @property const {
        return start_offset + bsize + 1;
    }

    public ushort bsize; /// total Block SIZE minus one

    public ushort cdata_size; /// compressed data size

    /// A buffer is used to reduce number of allocations.
    ///
    /// Its size is max(cdata_size, input_size)
    /// Initially, it contains compressed data, but is rewritten
    /// during decompressBgzfBlock -- indeed, who cares about
    /// compressed data after it has been uncompressed?
    public ubyte[] _buffer = void;

    /// If block has been already decompressed, result is undefined.
    public inout(ubyte[]) compressed_data() @property inout pure @safe nothrow {
        return _buffer[0 .. cast(size_t)cdata_size];
    }

    public uint crc32;
    public uint input_size; /// size of uncompressed data

    bool dirty;

    hash_t toHash() const pure @safe nothrow {
        assert(!dirty);
        return crc32;
    }

    bool opEquals(const ref BgzfBlock other) pure @safe nothrow {
        assert(!dirty);
        return opCmp(other) == 0;
    }

    int opCmp(const ref BgzfBlock other) const pure @safe nothrow {
        assert(!dirty);
        if (cdata_size < other.cdata_size)
            return -1;
        if (cdata_size > other.cdata_size)
            return 1;
        return std.algorithm.cmp(compressed_data, other.compressed_data);
    }
}

import std.stdio;

/**
  Struct representing decompressed BgzfBlock

  Start offset is needed to be able to tell current virtual offset,
  and yet be able to decompress blocks in parallel.
 */
struct DecompressedBgzfBlock {
  /* For the class version:
  this(ulong start, ulong end, ubyte[] buf) {
    start_offset = start;
    end_offset = end;
    decompressed_data = buf;
  }
  ~this() {
    stderr.writeln("destroy DecompressedBgzfBlock ",start_offset,":",end_offset," ",decompressed_data.sizeof);
  };
  */

  ulong start_offset;
  ulong end_offset;
  ubyte[] decompressed_data;
}

///
// alias Cache!(BgzfBlock, DecompressedBgzfBlock) BgzfBlockCache;

/// Function for BGZF block decompression.
/// Reuses buffer allocated for storing compressed data,
/// i.e. after execution buffer of the passed $(D block)
/// is overwritten with uncompressed data.
DecompressedBgzfBlock decompressBgzfBlock(BgzfBlock block)
{
    if (block.input_size == 0) {
      return DecompressedBgzfBlock(block.start_offset,
                                   block.start_offset + block.bsize + 1,
                                   cast(ubyte[])[]); // EOF marker
      // TODO: add check for correctness of EOF marker
    }

    /*
    if (cache !is null) {
        auto ptr = cache.lookup(block);
        if (ptr !is null)
            return *ptr;
    }
    */

    int err = void;

    // allocate buffer on the stack
    ubyte[BGZF_MAX_BLOCK_SIZE] uncompressed_buf = void;

    // check that block follows BAM specification
    enforce(block.input_size <= BGZF_MAX_BLOCK_SIZE,
            "Uncompressed block size must be within " ~
            to!string(BGZF_MAX_BLOCK_SIZE) ~ " bytes");

    // for convenience, provide a slice
    auto uncompressed = uncompressed_buf[0 .. block.input_size];

    // set input data
    bio.core.utils.zlib.z_stream zs;
    zs.next_in = cast(typeof(zs.next_in))block.compressed_data;
    zs.avail_in = to!uint(block.compressed_data.length);

    err = bio.core.utils.zlib.inflateInit2(&zs, /* winbits = */-15);
    if (err)
    {
        throw new ZlibException(err);
    }

    // uncompress it into a buffer on the stack
    zs.next_out = cast(typeof(zs.next_out))uncompressed_buf.ptr;
    zs.avail_out = block.input_size;

    err = bio.core.utils.zlib.inflate(&zs, Z_FINISH);
    switch (err)
    {
        case Z_STREAM_END:
            assert(zs.total_out == block.input_size);
            err = bio.core.utils.zlib.inflateEnd(&zs);
            if (err != Z_OK) {
                throw new ZlibException(err);
            }
            break;
        default:
            bio.core.utils.zlib.inflateEnd(&zs);
            throw new ZlibException(err);
    }

    assert(block.crc32 == crc32(0, uncompressed[]));

    /*
    if (cache !is null) {
        BgzfBlock compressed_bgzf_block = block;
        compressed_bgzf_block._buffer = block._buffer.dup;
        DecompressedBgzfBlock decompressed_bgzf_block;
        with (decompressed_bgzf_block) {
            start_offset = block.start_offset;
            end_offset = block.end_offset;
            decompressed_data = uncompressed[].dup;
        }
        cache.put(compressed_bgzf_block, decompressed_bgzf_block);
    }
    */

    // Now copy back to block._buffer, overwriting existing data.
    // It should have enough bytes already allocated.
    assert(block._buffer.length >= block.input_size);
    version(extraVerbose) {
        import std.stdio;
        stderr.writeln("[uncompressed] [write] range: ", block._buffer.ptr,
                       " - ", block._buffer.ptr + block.input_size);
    }
    block._buffer[0 .. block.input_size] = uncompressed[];
    block.dirty = true;

    auto decompressed = DecompressedBgzfBlock(block.start_offset, block.end_offset, block._buffer[0 .. block.input_size]);
    return decompressed;
}
