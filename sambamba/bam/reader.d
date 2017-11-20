/**
   New style BAM reader
*/

import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.exception;
import std.file;
import std.stdio;
import std.string;
import std.bitmanip;

struct Read {
}


struct Reader {

  immutable string filen;
  File f;

  this(string fn) {
    enforce(fn.isFile);
    filen = fn;

    ubyte[] tmp = new ubyte[4];
    f = File(fn,"r");
    auto magic = f.rawRead(tmp);
    assert(magic.length != 0);
    // note the magic number is encoded
    enforce(magic[0..4] == [0x1f,0x8b,0x08,0x04], "Invalid file format: expected file identifier BAM\\1 in "~filen);
  };

  int opApply(int delegate(ref int) operations) const  {
    return 0;
  }

}
