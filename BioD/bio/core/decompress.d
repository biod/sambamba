/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.decompress;

/**
   Streaming line reader which can be used for gzipped files. Note the
   current edition (still) uses the garbage collector. It may help to
   switch it off or to use the BioD decompressor used by bgzf.

   For a comparison with gzip a 2GB file decompressed with

   real    0m53.701s
   user    0m53.820s
   sys     0m0.572s

   while gzip took

   real    0m11.528s
   user    0m10.288s
   sys     0m0.936s

   So, that is something to aim for.

   Conversion can happen between different encodings, provided the
   line terminator is ubyte = '\n'. GzipbyLine logic is modeled on
   ByLineImpl and readln function from std.stdio.
*/

import std.algorithm;
// import std.concurrency;
import std.conv;
import std.exception;
import std.file;
import std.parallelism;
import std.stdio: File;
import std.zlib: UnCompress;

struct GzipbyLine(R) {

  File f;
  UnCompress decompress;
  R line;
  uint _bufsize;

  this(string gzipfn, uint bufsize=0x4000) {
    enforce(gzipfn.isFile);
    f = File(gzipfn,"r");
    decompress = new UnCompress();
    _bufsize = bufsize;
  }

  @disable this(this); // disable copy semantics;

  int opApply(scope int delegate(int line, R) dg) {

    int line = 0;
    // chunk_byLine takes a buffer and splits on \n.
    R chunk_byLine(R head, R rest) {
      auto split = findSplitAfter(rest,"\n");
      // If a new line is found split the in left and right.
      auto left = split[0]; // includes eol splitter
      auto right = split[1];
      if (left.length > 0) { // we have a match!
        dg(line++, head ~ left);
        return chunk_byLine([], right);
      }
      // no match
      return head ~ right;
    }

    R tail; // tail of previous buffer
    foreach (ubyte[] buffer; f.byChunk(_bufsize))
    {
      auto buf = cast(R)decompress.uncompress(buffer);
      tail = chunk_byLine(tail,buf);
    }
    if (tail.length > 0) dg(line++, tail);
    return 0;
  }
}


unittest {

  import std.algorithm.comparison : equal;

  // writeln("Testing GzipbyLine");
  int[] a = [ 1, 2, 4, 7, 7, 2, 4, 7, 3, 5];
  auto b = findSplitAfter(a, [7]);
  assert(equal(b[0],[1, 2, 4, 7]));
  assert(equal(b[1],[7, 2, 4, 7, 3, 5]));
  auto b1 = findSplitAfter(b[1], [7]);
  assert(equal(b1[0],[7]));
  assert(equal(b1[1],[2, 4, 7, 3, 5]));
  auto b2 = findSplitAfter([2, 4, 3], [7]);
  assert(equal(b2[0],cast(ubyte[])[]));
  assert(equal(b2[1],[2,4,3]));

  uint chars = 0;
  int lines = 0;
  /*
  foreach(line, ubyte[] s; GzipbyLine!(ubyte[])("test/data/BXD_geno.txt.gz")) {
    // test file contains 7320 lines 4707218 characters
    // write(cast(string)s);
    chars += s.length;
    lines = line;
  }
  */
  // These fail on recent versions of ldc
  // assert(lines == 7319,"genotype lines " ~ to!string(lines+1)); // fails with ldc2 < 1.10!
  // assert(chars == 4707218,"chars " ~ to!string(chars));
}

/**
   Mmfile threaded version of streaming line reader which can be used
   for gzipped files. Note the current edition is slower than
   GzipbyLine above and (still) uses the garbage collector. It may
   help to switch it off or to use the BioD decompressor used by bgzf.

   Conversion can happen between different encodings, provided the
   line terminator is ubyte = '\n'. GzipbyLine logic is modeled on
   ByLineImpl and readln function from std.stdio.
*/

import std.mmfile;
import core.thread;

struct GzipbyLineThreaded(R) {

  string fn;
  UnCompress decompress;
  R line;
  // Nullable!ubyte[] uncompressed_buf;
  uint _bufsize;

  this(string gzipfn, uint bufsize=0x4000) {
    enforce(gzipfn.isFile);
    fn = gzipfn;
    decompress = new UnCompress();
    _bufsize = bufsize;
  }

  @disable this(this); // disable copy semantics;

  int opApply(scope int delegate(int line, R) dg) {

    int line = 0;
    // chunk_byLine takes a buffer and splits on \n.
    R chunk_byLine(R head, R rest) {
      auto split = findSplitAfter(rest,"\n");
      // If a new line is found split the in left and right.
      auto left = split[0]; // includes eol splitter
      auto right = split[1];
      if (left.length > 0) { // we have a match!
        dg(line++, head ~ left);
        return chunk_byLine([], right);
      }
      // no match
      return head ~ right;
    }

    R decompressor(ubyte[] buffer) {
      return cast(R)decompress.uncompress(buffer);
    }

    auto mmf = new MmFile(fn);
    immutable mmf_length = mmf.length();
    long rest = mmf_length;
    R tail; // tail of previous buffer

    // Decompress the first chunk
    auto buffer1 = cast(ubyte[])mmf[0.._bufsize];
    rest -= buffer1.length;
    auto buf = decompressor(buffer1);

    uint chunknum = 1;
    while(rest>0) {
      // Get the next chunk
      ulong pos2 = (chunknum+1)*_bufsize;
      if (pos2 > mmf_length) pos2 = cast(ulong)mmf_length;
      auto buffer2 = cast(ubyte[])mmf[chunknum*_bufsize..mmf_length];
      rest -= buffer2.length;
      // Set up decompressing the next chunk
      auto t = task(&decompressor, buffer2);
      // auto t = task!decompressor(buffer2);
      t.executeInNewThread();
      // now invoke the delegate
      tail = chunk_byLine(tail,buf);
      buf = t.yieldForce();
      chunknum += 1;
    }
    tail = chunk_byLine(tail,buf);
    if (tail.length > 0) dg(line++, tail);
    return 0;
  }
}

unittest {
  int lines = 0;
  uint chars = 0;
  /*
  foreach(line, ubyte[] s; GzipbyLineThreaded!(ubyte[])("test/data/BXD_geno.txt.gz")) {
    // test file contains 7320 lines 4707218 characters
    // write(cast(string)s);
    chars += s.length;
    lines = line;
  }
  */
  /*
  These fail on recent versions of ldc
  assert(lines == 7319,"genotype lines " ~ to!string(lines+1));
  assert(chars == 4707218,"chars " ~ to!string(chars));
  */
}
