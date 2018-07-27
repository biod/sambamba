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
module sambamba.validate;

/**

   Validation.

   The new version is a prototype for new sambamba architecture using
   canonical D language features, including immutable and improved
   laziness and a more functional programming style. It should provide
   improved performance and minimize RAM use, as well as better
   composability. Also we are preparing it for CRAM input.

   Authors: Pjotr Prins

 */

import std.getopt;
import std.parallelism;
import std.range;
import std.stdio;
import std.typecons;

import bio2.bam.reader;
import bio2.bgzf;

void printUsage() {
  writeln("
Usage: sambamba validate [options] <input.bam> [<input2.bam> [...]]

       Validates a bam file.

Options:
         (...)
");
}

void info(string msg) {
  stderr.writeln("INFO: ",msg);
}

int validate_main(string[] args) {

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

    foreach (ref Nullable!ReadBlob read; BamReadBlobs(fn)) {
      auto pread = ProcessReadBlob(read);
      stdout.writeln(pread.toString);
    }
  }
  */
  return 0;
}
