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
module bio.core.bgzf.compress;

import bio.std.hts.bam.constants;
import bio.core.utils.zlib;

import std.array;
import std.system;
import std.bitmanip: nativeToLittleEndian;

/// Returns BGZF block containing compressed $(D chunk).
/// If $(D buffer) is provided, it will be used for storing the block.
///
/// Params:
///         chunk =  chunk of memory to be compressed
///         level =  compression level, see zlib documentation (https://github.com/madler/zlib/blob/master/zlib.h)
///         buffer = optional buffer which will be used for storing
///                  decompressed data
///
/// For uncompressed BAM output, use level = 0.
ubyte[] bgzfCompress(ubyte[] chunk, int level=-1, ubyte[] buffer=null)
in
{
    assert(-1 <= level && level <= 9);
}
body
{
    assert(bio.core.utils.zlib.compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);

    if (buffer is null) {
        buffer = uninitializedArray!(ubyte[])(BGZF_MAX_BLOCK_SIZE);
    } else {
        buffer.length = BGZF_MAX_BLOCK_SIZE;
    }

    // write header
    buffer[0 .. BLOCK_HEADER_LENGTH - ushort.sizeof] = BLOCK_HEADER_START[];

    bio.core.utils.zlib.z_stream zs;

    zs.zalloc = null;
    zs.zfree = null;

    zs.next_in  = cast(ubyte*)chunk.ptr;
    zs.avail_in = cast(uint)chunk.length;

    zs.next_out = buffer.ptr + BLOCK_HEADER_LENGTH;
    zs.avail_out = cast(int)(buffer.length - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH);

    auto err = bio.core.utils.zlib.deflateInit2(&zs, /* compression level */ level,
                                            /* deflated compression method */ Z_DEFLATED,
                                            /* winbits (no header) */ -15,
                                            /* memory usage level (default) */ 8,
                                            /* default compression strategy */ Z_DEFAULT_STRATEGY);
    if (err != Z_OK) {
        throw new ZlibException("deflateInit2", err);
    }

    err = bio.core.utils.zlib.deflate(&zs, Z_FINISH);
    if (err != Z_STREAM_END) {
        throw new ZlibException("deflate", err);
    }

    err = bio.core.utils.zlib.deflateEnd(&zs);
    if (err != Z_OK) {
        throw new ZlibException("deflateEnd", err);
    }

    // almost done, update buffer length
    buffer.length = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

    // Write (block length - 1) in BC subfield.
    // Why -1? To fit the value into 2 bytes (it's assumed to be in range 1-65536).
    ushort len = cast(ushort)(buffer.length - 1);
    buffer[BLOCK_HEADER_LENGTH - 2 .. $][0 .. 2] = nativeToLittleEndian(len);

    // Write the footer
    buffer[$ - 8 .. $ - 4] = nativeToLittleEndian(crc32(0, chunk));
    buffer[$- 4 .. $] = nativeToLittleEndian(cast(uint)chunk.length);
    return buffer;
}
