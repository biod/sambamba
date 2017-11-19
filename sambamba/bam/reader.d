/**
   New style BAM reader
*/

import std.exception;
import std.file;

struct Read {
}

struct Reader {

  immutable string _filename;

  this(string fn) {
    enforce(fn.isFile);
    _filename = fn;
  };

  int opApply(int delegate(ref int) operations) const  {
    return 0;
  }

}
