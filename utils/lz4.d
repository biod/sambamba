/*
    This file is part of Sambamba.
    Copyright (C) 2012-2015    Artem Tarasov <lomereiter@gmail.com>

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
module utils.lz4;
import std.stdio, undead.stream, std.conv;

extern(C) {
  alias size_t LZ4F_errorCode_t;
  uint LZ4F_isError(LZ4F_errorCode_t code);

  immutable(char)* LZ4F_getErrorName(LZ4F_errorCode_t code); /* return error code string; useful for debugging */

  enum blockSizeID_t { LZ4F_default=0, max64KB=4, max256KB=5, max1MB=6, max4MB=7 }
  enum blockMode_t { blockLinked=0, blockIndependent}
  enum contentChecksum_t { noContentChecksum=0, contentChecksumEnabled }
  enum frameType_t { LZ4F_frame=0, skippableFrame }

  struct LZ4F_frameInfo_t {
    blockSizeID_t blockSizeID; /* max64KB, max256KB, max1MB, max4MB ; 0 == default */
    blockMode_t blockMode; /* blockLinked, blockIndependent ; 0 == default */
    contentChecksum_t contentChecksumFlag; /* noContentChecksum, contentChecksumEnabled ; 0 == default */
    frameType_t frameType; /* LZ4F_frame, skippableFrame ; 0 == default */
    ulong contentSize; /* Size of uncompressed (original) content ; 0 == unknown */
    uint[2] reserved; /* must be zero for forward compatibility */
  }

  struct LZ4F_preferences_t {
    LZ4F_frameInfo_t frameInfo;
    uint compressionLevel; /* 0 == default (fast mode); values above 16 count as 16 */
    uint autoFlush; /* 1 == always flush (reduce need for tmp buffer) */
    uint[4] reserved; /* must be zero for forward compatibility */
  }

  size_t LZ4F_compressFrameBound(size_t srcSize, const LZ4F_preferences_t* preferencesPtr);
  size_t LZ4F_compressFrame(void* dstBuffer, size_t dstMaxSize, const void* srcBuffer, size_t srcSize, const LZ4F_preferences_t* preferencesPtr);

  alias void* LZ4F_compressionContext_t;
  struct LZ4F_compressOptions_t {
    uint stableSrc; /* 1 == src content will remain available on future calls to LZ4F_compress(); avoid saving src content within tmp buffer as future dictionary */
    uint[3] reserved;
  }

  enum LZ4F_VERSION = 100;
  LZ4F_errorCode_t LZ4F_createCompressionContext(LZ4F_compressionContext_t* cctxPtr, uint version_);
  LZ4F_errorCode_t LZ4F_freeCompressionContext(LZ4F_compressionContext_t cctx);
  size_t LZ4F_compressBegin(LZ4F_compressionContext_t cctx, void* dstBuffer, size_t dstMaxSize, const LZ4F_preferences_t* prefsPtr);
  size_t LZ4F_compressBound(size_t srcSize, const LZ4F_preferences_t* prefsPtr);
  size_t LZ4F_compressUpdate(LZ4F_compressionContext_t cctx, void* dstBuffer, size_t dstMaxSize, const void* srcBuffer, size_t srcSize, const LZ4F_compressOptions_t* cOptPtr);
  size_t LZ4F_flush(LZ4F_compressionContext_t cctx, void* dstBuffer, size_t dstMaxSize, const LZ4F_compressOptions_t* cOptPtr);
  size_t LZ4F_compressEnd(LZ4F_compressionContext_t cctx, void* dstBuffer, size_t dstMaxSize, const LZ4F_compressOptions_t* cOptPtr);
  alias void* LZ4F_decompressionContext_t;
  struct LZ4F_decompressOptions_t {
    uint stableDst; /* guarantee that decompressed data will still be there on next function calls (avoid storage into tmp buffers) */
    uint[3] reserved;
  }
  LZ4F_errorCode_t LZ4F_createDecompressionContext(LZ4F_decompressionContext_t* dctxPtr, uint version_);
  LZ4F_errorCode_t LZ4F_freeDecompressionContext(LZ4F_decompressionContext_t dctx);
  size_t LZ4F_getFrameInfo(LZ4F_decompressionContext_t dctx, LZ4F_frameInfo_t* frameInfoPtr, const void* srcBuffer, size_t* srcSizePtr);
  size_t LZ4F_decompress(LZ4F_decompressionContext_t dctx, void* dstBuffer, size_t* dstSizePtr, const void* srcBuffer, size_t* srcSizePtr, const LZ4F_decompressOptions_t* dOptPtr);

}

class LZ4Exception : Exception {
  this(string description, LZ4F_errorCode_t errorCode) {
    super(description ~ " : " ~ LZ4F_getErrorName(errorCode).to!string);
  }

  this(string description) {
    super(description);
  }
}

class LZ4Compressor {
  private {
    LZ4F_preferences_t prefs;
    ubyte[] in_buff;
    ubyte[] out_buff;
    LZ4F_compressionContext_t ctx;

    LZ4F_compressionContext_t createCompressionContext() {
      LZ4F_compressionContext_t ctx;
      auto code = LZ4F_createCompressionContext(&ctx, LZ4F_VERSION);
      if (LZ4F_isError(code))
        throw new LZ4Exception("Failure in LZ4F_createCompressionContext", code);
      return ctx;
    }

    void freeCompressionContext(LZ4F_compressionContext_t ctx) {
      auto code = LZ4F_freeCompressionContext(ctx);
      if (LZ4F_isError(code))
        throw new LZ4Exception("Failed to free LZ4F compression context", code);
    }
  }

  this() {
    prefs.autoFlush = 1;
    prefs.frameInfo.blockMode = blockMode_t.blockIndependent;
    prefs.frameInfo.blockSizeID = blockSizeID_t.max256KB;
    prefs.frameInfo.contentChecksumFlag = contentChecksum_t.contentChecksumEnabled;

    int block_size = (1 << (8 + (2 * prefs.frameInfo.blockSizeID)));
    in_buff.length = block_size;
    out_buff.length = LZ4F_compressBound(block_size, &prefs);
  }

  private {
    size_t compressBegin() {
      auto sz = LZ4F_compressBegin(ctx, out_buff.ptr, out_buff.length, &prefs);
      if (LZ4F_isError(sz))
        throw new LZ4Exception("Failure in LZ4F_compressBegin", sz);
      return sz;
    }

    size_t compressUpdate(ubyte[] block) {
      auto sz = LZ4F_compressUpdate(ctx, out_buff.ptr, out_buff.length, block.ptr, block.length, null);
      if (LZ4F_isError(sz))
        throw new LZ4Exception("Failure in LZ4F_compressUpdate", sz);
      return sz;
    }

    size_t compressEnd() {
      auto sz = LZ4F_compressEnd(ctx, out_buff.ptr, out_buff.length, null);
      if (LZ4F_isError(sz))
        throw new LZ4Exception("Failure in LZ4F_compressEnd", sz);
      return sz;
    }
  }

  void compress(ubyte[] function(ubyte[], void*) read_block, void* data,
                std.stdio.File output_file,
                int compressionLevel=0)
  {
    prefs.compressionLevel = compressionLevel;

    ctx = createCompressionContext();
    assert(ctx !is null);

    auto sz = compressBegin();
    output_file.rawWrite(out_buff[0 .. sz]);

    while (true) {
      ubyte[] block = read_block(in_buff, data);
      if (block.length == 0) break;
      sz = compressUpdate(block);
      output_file.rawWrite(out_buff[0 .. sz]);
    }

    sz = compressEnd();
    output_file.rawWrite(out_buff[0 .. sz]);

    freeCompressionContext(ctx);
  }

  private {
    static ubyte[] read_block_from_file(ubyte[] buf, void* data) {
      auto input_file = cast(std.stdio.File*)data;
      return input_file.rawRead(buf);
    }

    import std.algorithm : min;
    static ubyte[] read_block_from_array(ubyte[] buf, void* data) {
      auto arr = cast(ubyte[]*)data;
      auto n = min(buf.length, (*arr).length);
      buf[0 .. n] = (*arr)[0 .. n];
      *arr = (*arr)[n .. $];
      return buf[0 .. n];
    }
  }
  void compress(std.stdio.File input_file, std.stdio.File output_file,
                int compressionLevel=0) {
    compress(&read_block_from_file, &input_file, output_file, compressionLevel);
  }

  void compress(ubyte[] data, std.stdio.File output_file,
                int compressionLevel=0) {
    compress(&read_block_from_array, &data, output_file, compressionLevel);
  }
}

private {
  LZ4F_decompressionContext_t createDecompressionContext() {
    LZ4F_decompressionContext_t ctx;
    auto code = LZ4F_createDecompressionContext(&ctx, LZ4F_VERSION);
    if (LZ4F_isError(code))
      throw new LZ4Exception("Failure in LZ4F_createDecompressionContext", code);
    return ctx;
  }

  void freeDecompressionContext(LZ4F_compressionContext_t ctx) {
    auto code = LZ4F_freeDecompressionContext(ctx);
    if (LZ4F_isError(code))
      throw new LZ4Exception("Failed to free LZ4F decompression context", code);
  }
}

struct LZ4File {
  private {
    std.stdio.File input_file;
    LZ4F_decompressionContext_t ctx;

    size_t bytes_read, bytes_written;

    void decompress(ubyte[] block, ubyte[] out_buff) {
        bytes_read = block.length;
        bytes_written = out_buff.length;
        auto code = LZ4F_decompress(ctx, out_buff.ptr, &bytes_written, block.ptr, &bytes_read, null);
        if (LZ4F_isError(code))
          throw new LZ4Exception("Failure in LZ4F_decompress", code);
    }

    ubyte[] in_buff;
    ubyte[] block;
  }

  this(string filename) {
    input_file = std.stdio.File(filename);

    in_buff.length = 256 << 10;

    ctx = createDecompressionContext();
    assert(ctx !is null);

    auto header = input_file.rawRead(in_buff[0 .. 4]);
    assert(header.length == 4);

    ubyte[256] tmp;
    decompress(in_buff[0 .. 4], tmp[]); // header
    assert(bytes_read == 4);
    assert(bytes_written == 0);
  }

  void close() {
    freeDecompressionContext(ctx);
    input_file.close();
  }

  ubyte[] rawRead(ubyte[] buffer) {
    while (true) {
      decompress(block, buffer);
      block = block[bytes_read .. $];
      if (bytes_written > 0)
        break;
      if (block.length == 0 && !readNextBlock())
        break;
    }
    return buffer[0 .. bytes_written];
  }

  private bool readNextBlock() {
      block = input_file.rawRead(in_buff);
      return block.length > 0;
  }
}

class LZ4Decompressor {
  private {
    ubyte[] in_buff;
    ubyte[] out_buff;
  }

  this() {
    in_buff.length = 256 << 10;
    out_buff.length = 256 << 10;
  }

  void decompress(undead.stream.InputStream input_file,
                  std.stdio.File output_file) {
    size_t bytes_read, bytes_written;

    auto ctx = createDecompressionContext();
    assert(ctx !is null);

    void decompress(ubyte[] block) {
      bytes_read = block.length;
      bytes_written = out_buff.length;
      auto code = LZ4F_decompress(ctx, out_buff.ptr, &bytes_written, block.ptr, &bytes_read, null);
      if (LZ4F_isError(code))
        throw new LZ4Exception("Failure in LZ4F_decompress", code);
    }

    input_file.readExact(in_buff.ptr, 4);

    decompress(in_buff[0 .. 4]); // header
    assert(bytes_read == 4);
    assert(bytes_written == 0);

    while (true) {
      size_t raw_read = input_file.read(in_buff);
      if (raw_read == 0) break;
      auto block = in_buff[0 .. raw_read];
      while (block.length > 0) {
        decompress(block);
        block = block[bytes_read .. $];
        if (bytes_written > 0)
          output_file.rawWrite(out_buff[0 .. bytes_written]);
      }
    }
    freeDecompressionContext(ctx);
  }
}

int lz4compress_main() {
  auto compressor = new LZ4Compressor();
  compressor.compress(stdin, stdout);
  return 0;
}
