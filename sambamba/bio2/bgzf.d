
import std.bitmanip;
import std.conv;
import std.exception;
import std.file;
import std.stdio;

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

  size_t get_compressed_block(size_t fpos) {
    void throwBgzfException(string msg) {
        throw new BgzfException("Error reading BGZF block starting in "~filen~" @ " ~
                                to!string(fpos) ~ ": " ~ msg);
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

    auto start_offset = fpos;
    write(".");
    f.seek(fpos);

    auto magic = f.rawRead(new ubyte[4]);
    enforce(magic[0..4] == BGZF_MAGIC, "Invalid file format: expected bgzf magic number in "~filen~" @ " ~ to!string(fpos));
    enforce(magic.length == 4,"Premature end of file in ~",filen);
    try {
      f.rawRead(new ubyte[uint.sizeof + 2 * ubyte.sizeof]); // skip gzip info
      ushort gzip_extra_length = read_ushort();
      writeln("gzip_extra_length=",gzip_extra_length);
      writeln([f.tell, fpos, gzip_extra_length]);
      immutable fpos1 = f.tell;
      ushort bsize = 0;
      while (f.tell < fpos1 + gzip_extra_length) {
        auto subfield_id1 = read_ubyte();
        auto subfield_id2 = read_ubyte();
        auto subfield_len = read_ushort();
        if (subfield_id1 == BAM_SI1 && subfield_id2 == BAM_SI2) {
          // BC identifier
          enforce(gzip_extra_length == 6);
          // FIXME: check there is only one BC block
          bsize = read_ushort(); // BLOCK size minus 1
          enforce(subfield_len == bsize.sizeof, "BC subfield len should be 2");
          writeln("bsize=",bsize);
          break;
        }
        else {
          f.seek(subfield_len,SEEK_CUR);
        }
      }
      if (bsize==0) throwBgzfException("block size not found");
      auto cdata_size = bsize - gzip_extra_length - 19;
      if (cdata_size > BGZF_MAX_BLOCK_SIZE) throwBgzfException("compressed size larger than allowed");
      // @@ HERE
    } catch (Exception e) { throwBgzfException("File error in "~filen~": " ~ e.msg); }
    return f.tell()-start_offset;
  }

  string blocks() {
    string ret = "yes";
    auto fpos = 0;
    while (!f.eof()) {
      auto block_length = get_compressed_block(fpos);
      fpos += block_length;
    }
    return ret;
  }
}
