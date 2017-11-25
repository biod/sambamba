
import std.bitmanip;
import std.conv;
import std.exception;
import std.file;
import std.stdio;
import std.typecons;
import std.zlib;

import bio.bam.constants;
import bio.core.bgzf.block;
import bio.core.bgzf.constants;

class BgzfException : Exception {
    this(string msg) { super(msg); }
}

struct BgzfReader {
  string filen;
  File f;

  this(string fn) {
    enforce(fn.isFile);
    filen = fn;
    f = File(fn,"r");
  }

  /**
   * Returns new file position. Zero when done
   */
  Nullable!size_t get_compressed_block(Nullable!size_t fpos) {
    void throwBgzfException(string msg, string file = __FILE__, int line = __LINE__) {
        throw new BgzfException("Error reading BGZF block starting in "~filen~" @ " ~
                                to!string(fpos) ~ " (" ~ file ~ ":" ~ to!string(line) ~ "): " ~ msg);
    }
    void enforce1(bool check, lazy string msg, string file = __FILE__, int line = __LINE__) {
      if (!check)
        throwBgzfException(msg,file,line);
    }
    ubyte read_ubyte() {
      ubyte[1] ubyte1; // read buffer
      immutable ubyte[1] buf = f.rawRead(ubyte1);
      return buf[0];
    }
    ushort read_ushort() {
      ubyte[2] ubyte2; // read buffer
      immutable ubyte[2] buf = f.rawRead(ubyte2);
      return littleEndianToNative!ushort(buf);
    }
    auto read_uint() {
      ubyte[4] ubyte4; // read buffer
      immutable ubyte[4] buf = f.rawRead(ubyte4);
      return littleEndianToNative!uint(buf);
    }

    immutable start_offset = fpos;
    write(".");
    f.seek(fpos);

    ubyte[4] ubyte4;
    auto magic = f.rawRead(ubyte4);
    enforce1(magic.length == 4, "Premature end of file");
    enforce1(magic[0..4] == BGZF_MAGIC,"Invalid file format: expected bgzf magic number");
    try {
      ubyte[uint.sizeof + 2 * ubyte.sizeof] skip;
      f.rawRead(skip); // skip gzip info
      ushort gzip_extra_length = read_ushort();
      // writeln("gzip_extra_length=",gzip_extra_length);
      immutable fpos1 = f.tell;
      size_t bsize = 0;
      while (f.tell < fpos1 + gzip_extra_length) {
        immutable subfield_id1 = read_ubyte();
        immutable subfield_id2 = read_ubyte();
        immutable subfield_len = read_ushort();
        if (subfield_id1 == BAM_SI1 && subfield_id2 == BAM_SI2) {
          // BC identifier
          enforce(gzip_extra_length == 6);
          // FIXME: always picks first BC block
          bsize = 1+read_ushort(); // BLOCK size
          enforce1(subfield_len == 2, "BC subfield len should be 2");
          // writeln("bsize=",bsize);
          break;
        }
        else {
          f.seek(subfield_len,SEEK_CUR);
        }
      }
      enforce1(bsize!=0,"block size not found");
      f.seek(fpos1+gzip_extra_length); // skip any extra subfields - note we don't check for second BC
      immutable compressed_size = bsize - 1 - gzip_extra_length - 19;
      enforce1(compressed_size <= BGZF_MAX_BLOCK_SIZE, "compressed size larger than allowed");

      stderr.writeln("[compressed] size ", compressed_size, " bytes starting block @ ", start_offset);
      ubyte[BGZF_MAX_BLOCK_SIZE] stack_buffer;
      auto buffer = stack_buffer[0..compressed_size];
      auto compressed_buf = f.rawRead(buffer);

      immutable crc32 = read_uint();
      immutable uncompressed_size = read_uint();
      stderr.writeln("[uncompressed] size ",uncompressed_size);

      if (uncompressed_size == 0) {
        // check for eof marker, rereading block header
        auto lastpos = f.tell();
        f.seek(start_offset);
        ubyte[28] buf;
        f.rawRead(buf);
        f.seek(lastpos);
        if (buf == [31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0]) return Nullable!size_t();
      }

      auto uncompressed_buf = uncompress(compressed_buf,compressed_size,-15);
      assert(crc32 == std.zlib.crc32(0, uncompressed_buf[]));

    } catch (Exception e) { throwBgzfException("File error in "~filen~": " ~ e.msg); }
    return Nullable!size_t(f.tell());
  }

  string blocks() {
    string ret = "yes";
    Nullable!size_t fpos = 0;
    while (!f.eof()) {
      auto new_fpos = get_compressed_block(fpos);
      if (new_fpos.isNull) {
        writeln(to!string(fpos) ~ "bytes");
        break;
      }
      fpos = new_fpos;
    }
    return ret;
  }
}
