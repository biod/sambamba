/*
    This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/
module sambamba.markdup2;

/**

   Markdup2 is a new version of sambamba markdup (under development).

   The new version is a prototype for new sambamba architecture using
   canonical D language features, including immutable and improved
   laziness and a more functional programming style. It should provide
   improved performance and minimize RAM use, as well as better
   composability. Also we are preparing it for CRAM input.

   The initial version is a large data markdup which was previously invoked as

     sambamba markdup --hash-table-size=2097152 --overflow-list-size=160000 --io-buffer-size=1024 in.bam out.bam

   Authors: Pjotr Prins

 */

import std.getopt;
import std.parallelism;
import std.range;
import std.stdio;

import bio.std.experimental.hts.bam.reader;
import bio.std.experimental.hts.bgzf;

void printUsage() {
  writeln("
Usage: sambamba markdup2 [options] <input.bam> [<input2.bam> [...]]

       By default, marks duplicates without removing them.
       Writes to stdout if no output file is given.

Options: -r, --remove-duplicates   remove duplicates instead of just marking them (nyi)
         -o, --output-filename fn  write to output file (bam format) (nyi)
");
}

void info(string msg) {
  stderr.writeln("INFO: ",msg);
}

int markdup_main(string[] args) {
  bool remove_duplicates;
  string outfn;
  getopt(args,
         std.getopt.config.caseSensitive,
         "remove-duplicates|r", &remove_duplicates,
         "output-filename|o", &outfn);

  if (args.length < 2) {
    printUsage();
    return 1;
  }

  info("Reading input files");

  auto infns = args[1..$];

  // auto taskpool = new TaskPool();
  // scope(exit) taskpool.stop();

  /*
  foreach (string fn; infns) {
    stderr.writeln(fn);
    foreach (ref ReadBlob read; BamReadBlobs(fn)) {
      auto pread = ProcessReadBlob(read); // FIXME we don't need ProcessRead here
      writeln(pread.toString, ",", pread.start_pos, ",", pread.end_pos);
      // stdout.write(read);
    }
  }
  */
  return 0;
}
