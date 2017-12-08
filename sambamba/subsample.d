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
module sambamba.subsample;

/**

   Subsampling.

   Subsampling limits the read depth to a certain threshold. This is
   increasingly important with large sequencing efforts where high
   depth can be reached, especially in non-informative regions - i.e.,
   with high repeats. A good subsample method does not have to be
   exact (we can vary around the maximum read depth) but has to be
   reproducible and limit effect on down stream variant calling.

   The first algorithm 'hash1' is the same one used in
   VariantBam::SubSampleWrite. It simply takes the coverages at the
   beginning and end of a read, takes the maximum and drops reads
   based on a Hash computation (you end up with an approximate number
   of reads around max_depth). This algorithm is reproducible but does
   not consider other factors and read pairs.

   Authors: Pjotr Prins and Brad Chapman

 */

import std.getopt;
import std.parallelism;
import std.range;
import std.stdio;

import sambamba.bio2.bam.reader;
import sambamba.bio2.bgzf;
import sambamba.bio2.constants;
import sambamba.bio2.pileup;

void printUsage() {
  writeln("
Usage: sambamba subsample [options] <input.bam> [<input2.bam> [...]]

       Subsample a bam file.

Options:

         --type [hash1]  Algorithm for subsampling (hash1, default is none)
");
}

void info(string msg) {
  stderr.writeln("INFO: ",msg);
}

/**
   For every read track position
*/
struct ReadInfo {
  RefId ref_id;
  GenomePos start, stop;

  this(ref ProcessRead2 read) {
    ref_id = read.ref_id;
    start = read.start_pos;
    stop = read.end_pos;
  }
}

int subsample_main(string[] args) {

  if (args.length < 2) {
    printUsage();
    return 1;
  }

  info("Reading input files");

  auto infns = args[1..$];

  // auto taskpool = new TaskPool();
  // scope(exit) taskpool.stop();

  auto pileup = new PileUp!ReadInfo();

  foreach (string fn; infns) {
    stderr.writeln(fn);

    foreach (ref Read2 read; BamReader2(fn)) {
      auto pread = ProcessRead2(read);
      auto r = ReadInfo(pread);
      pileup.push(r);
      writeln(pread.toString, ",", pread.start_pos, ",", pread.end_pos);
    }
  }
  return 0;
}
