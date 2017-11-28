
import std.bitmanip;
import std.conv;
import std.exception;
import std.file;
import std.stdio;
import std.typecons;
import std.zlib : calc_crc32 = crc32, ZlibException;

import bio.bam.constants;
import bio.core.bgzf.block;
import bio.core.bgzf.constants;
import bio.core.utils.zlib : inflateInit2, inflate, inflateEnd, Z_OK, Z_FINISH, Z_STREAM_END;

class BgzfException : Exception {
    this(string msg) { super(msg); }
}

alias Nullable!size_t FilePos;
alias immutable(uint) CRC32;

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
  assert(uncompressed_buf.length == BGZF_MAX_BLOCK_SIZE);
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

  return cast(ubyte[])uncompressed_buf;
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
    enforce1(compressed_size <= BGZF_MAX_BLOCK_SIZE, "compressed size larger than allowed");

    stderr.writeln("[compressed] size ", compressed_size, " bytes starting block @ ", report_fpos);
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
    stderr.writeln("[uncompressed] size ",uncompressed_size);
    return tuple(compressed_buf,uncompressed_size,crc32);
  }

  /**
   * Returns new tuple of the new file position, the compressed buffer and
   * the CRC32 o the uncompressed data. file pos is NULL when done
   */
  Tuple!(FilePos,ubyte[],size_t,CRC32) read_compressed_block(FilePos fpos, ubyte[] buffer) {
    immutable start_offset = fpos;
    try {
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
        if (buf == BGZF_EOF) return tuple(FilePos(),compressed_buf,cast(ulong)0,crc32);
      }
      return tuple(FilePos(f.tell()),compressed_buf,cast(size_t)uncompressed_size,crc32);
    } catch (Exception e) { throwBgzfException(e.msg,e.file,e.line); }
    assert(0); // never reached
  }
}

struct BgzfBlocks {

  BgzfReader bgzf;

  this(string fn) {
    bgzf = BgzfReader(fn);
  }

  int opApply(scope int delegate(ubyte[]) dg) {
    FilePos fpos = 0;

    try {
      while (!fpos.isNull) {
        ubyte[BGZF_MAX_BLOCK_SIZE] stack_buffer;
        auto res = bgzf.read_compressed_block(fpos,stack_buffer);
        fpos = res[0]; // point fpos to next block
        if (fpos.isNull) break;

        auto compressed_buf = res[1]; // same as stack_buffer
        auto uncompressed_size = res[2];
        auto crc32 = res[3];
        ubyte[BGZF_MAX_BLOCK_SIZE] uncompressed_buf;
        // call delegated function with new block
        dg(deflate(uncompressed_buf,compressed_buf,uncompressed_size,crc32));
      }
    } catch (Exception e) { bgzf.throwBgzfException(e.msg,e.file,e.line); }
    return 0;
  }
}
