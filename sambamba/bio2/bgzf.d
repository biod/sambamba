
import std.conv;
import std.exception;
import std.file;
import std.stdio;

import bio.core.bgzf.block;

struct BgzfReader {
  string filen;
  File f;

  this(string fn) {
    enforce(fn.isFile);
    filen = fn;
    f = File(fn,"r");
  }

  size_t get_compressed_block(size_t fpos) {
    auto fpos1 = fpos;
    write(".");
    f.seek(fpos);
    auto magic = f.rawRead(new ubyte[4]);
    enforce(magic[0..4] == BGZF_MAGIC, "Invalid file format: expected bgzf magic number in "~filen~" @ " ~ to!string(fpos));
    enforce(magic.length == 4,"Premature end of file in ~",filen);
    return f.tell()-fpos1;
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
