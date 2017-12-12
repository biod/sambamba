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

import std.experimental.logger;
import std.getopt;
import std.parallelism;
import std.range;
import std.stdio;
import std.typecons;

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
         --logging type  Set logging to debug|info|warning|critical
");
}

/**
   For every read track position
*/
struct ReadInfo {
  RefId ref_id;
  GenomePos start_pos, end_pos;

  this(ref ProcessReadBlob read) {
    ref_id = read.ref_id;
    start_pos = read.start_pos;
    end_pos = read.end_pos;
  }
}

int subsample_main2(string[] args) {
  globalLogLevel(LogLevel.trace); // debug level

  if (args.length < 2) {
    printUsage();
    return 1;
  }

  auto infns = args[1..$];

  auto pileup = new PileUp!ReadInfo();
  Nullable!ReadInfo prev; // keep track of previous reads

  foreach (string fn; infns) {
    stderr.writeln(fn);

    foreach (ref ReadBlob read; BamReadBlobStream(fn)) {
      auto pread = ProcessReadBlob(read); // FIXME we don't need ProcessRead here
      // Read ahead until the window is full (FIXME)
      auto r = ReadInfo(pread);
      pileup.push(r);
      writeln(pread.toString, ",", pread.start_pos, ",", pread.end_pos);
      if (!prev.isNull) {
        // Remove reads that have gone out of the window (FIXME)
        /*
        pileup.delete_if( stacked_read =>
                          stacked_read.ref_id != r.ref_id || stacked_read.end_pos < r.begin_pos
                          );
        */
      }
      // prev = r;
    }
  }
  return 0;
}


int subsample_main(string[] args) {
  globalLogLevel(LogLevel.trace); // debug level

  if (args.length < 2) {
    printUsage();
    return 1;
  }

  auto infns = args[1..$];

  auto pileup = new PileUp!ProcessReadBlob(1000);

  foreach (string fn; infns) {
    stderr.writeln(fn);

    auto stream = BamReadBlobStream(fn);

    int popped = 0;
    // get the first two reads
    auto lread = ProcessReadBlob(stream.read);
    auto lread_idx = pileup.push(lread);

    while (!stream.empty) {
      lread = pileup.read_at(lread_idx);
      auto rread = lread;
      // Fill ring buffer ahead until the window is full
      while (!stream.empty && lread.ref_id == rread.ref_id && rread.start_pos < lread.end_pos+1) {
        rread = ProcessReadBlob(stream.read);
        rread_idx = pileup.push(rread);
      }
      // Now we have a pileup and we can check this read
      writeln("---> stopped pileup at ",lread.ref_id," ",lread.start_pos,":",lread.end_pos," ",rread.start_pos,":",rread.end_pos);

      ulong depth = pileup.ldepth(lread.start_pos,lread_idx);
      writeln("*** ",lread.start_pos," ",lread_idx," Read depth is ",depth);
      // Remove the current read
      pileup.popFront();
      popped++;
      writeln("Ring buffer size ",pileup.ring.length);
      writeln("popped ",popped);
      writeln("pushed ",pileup.ring.pushed," popped ",pileup.ring.popped);

      lread_idx += 1;

      // if (!prev.isNull) {
        // Remove reads that have gone out of the window (FIXME)
        /*
        pileup.delete_if( stacked_read =>
                          stacked_read.ref_id != r.ref_id || stacked_read.end_pos < r.begin_pos
                          );
        */
      // }
      // prev = r;
    }
  }
  return 0;
}
