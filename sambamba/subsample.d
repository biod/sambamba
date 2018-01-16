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

import core.memory : GC;
import std.algorithm.comparison : max;
import std.conv;
import std.experimental.logger;
import std.exception;
import std.getopt;
import std.parallelism;
import std.range;
import std.stdio;
import std.typecons;

import sambamba.bio2.bam.reader;
import sambamba.bio2.bam.writer;
import sambamba.bio2.bgzf;
import sambamba.bio2.hashing;
import sambamba.bio2.constants;
import sambamba.bio2.pileup;
import sambamba.bio2.reads;

import bio.core.utils.exception;

void printUsage() {
  writeln("
Usage: sambamba subsample [options] <input.bam> [<input2.bam> [...]]

       Subsample a bam file.

Options:

         --type [fasthash]   Algorithm for subsampling (fasthash, default is none)
         --max-cov [depth]   Maximum coverage (approx)
         -o, --output fn     Set output file (default stdout)
         -r, --remove        Remove over sampled reads from output

         --logging type   Set logging to debug|info|warning|critical -nyi

Examples:

       sambamba subsample --type=fasthash input.bam -ooutput.bam
");
}

enum RState {
  unknown,     // unknown/unmapped read (always write)
  keep,        // keep in set (allways write)
  drop,        // drop from set (no write or mark as bad quality)
  dirty        // delete from ring buffer (done)
}

// ReadState keeps track of the state of a processed Read. This state
// is maintained on the ringbuffer. We may change this design later.

struct ReadState {
  ProcessReadBlob read;
  RState state;

  this(ProcessReadBlob _r) {
    enforce(!_r.isNull);
    read = _r;
    state = RState.unknown;
  }

  // @disable this(this); // disable copy semantics;

  @property cleanup() {
    assert(is_dirty);
    read.cleanup;
  }
  @property ref ProcessReadBlob get() {
    return read;
  }
  @property void set_keep() {
    state = RState.keep;
  }
  @property void set_drop() {
    state = RState.drop;
  }
  @property void set_dirty() {
    state = RState.dirty;
  }
  @property bool is_dropped() {
    return state == RState.drop;
  }
  @property bool is_dirty() {
    return state == RState.dirty;
  }
}


/**
   Implementation of fasthash which is the simplest implementation,
   comparable to that of others.

   While reads stream in they get piled up in a ringbuffer. The
   ringbuffer gets filled ahead until the read that is no longer
   overlapping, creating a window from leftmost to rightmost:

                                                          r-------------
                                               ----y---------------
                                        -----------y--------
                                 x=================y
                            -----x-------------
                   --------------x----
              l------------------x-----
          leftmost             start_pos       end_pos  rightmost

   once the reads have been counted it moves to the next read, reads
   ahead and destroys the leftmost reads that have gone out of the
   window.

   Some read stacks are (theoretically) unlimited in size. Therefore
   we stop processing them at (say) 20x the max_cov (max_reached). At
   that point the reader goes into a separate mode, continuing to read
   and write on the fly without filling the ringbuffer.

   Depth is cached in a separate ringbuffer at start positions. The
   depth at the end_pos is inferred/estimated from this.

   In pseudo code:

   each read
     read+write while is_ignore == unmapped/bad quality/duplicates
     read ahead while overlapping
       if max_coverage reached
         go into separate blocking mode
         write out stack
     compute depth
     store depth in cache
     mark read for keep/drop
     destroy reads that are no longer overlapping (reaper)

   Pretty simple. In the future we may add an extra layer for
   processing, but I don't think it will help much.

*/

// ---- Move through reads that are ignored. Modifies pileup
void foreach_invalid_read(ref BamReadBlobStream reader, void delegate(ProcessReadBlob) dg) {
  auto read = ProcessReadBlob(reader.read);
  // auto first = pileup.push(ReadState(first_read)); // the current pointer maintains state with the pileup
  while(!read.isNull && read.is_unmapped) {
    dg(read);
    read = ProcessReadBlob(reader.read);
  }
}

  // auto read = pileup.get(current);
  /*
    assert(!current.isNull);
    while (current.is_unmapped2) {
    // we hit an unmapped set, need to purge (this won't work on threads)
    // writeln("Skip unmapped read");
    reap();
    current = ProcessReadBlob(stream.read);
    if (current.isNull)
    break;
    current_idx = pileup.push(ReadState(current));
    current = pileup.read_at(current_idx).get;
    // rightmost = current;
    rightmost_idx = current_idx;
    }
    assert(current.is_mapped2);
  */

int subsample_main(string[] args) {
  bool remove = false;
  globalLogLevel(LogLevel.trace); // debug level
  int max_cov = 0;

  if (args.length < 2) {
    printUsage();
    return 1;
  }

  string outputfn;
  string type;

  getopt(args,
         std.getopt.config.caseSensitive,
         "type", &type,
         "max-cov", &max_cov,
         "output|o", &outputfn,
         "remove|r", &remove,
         );

  enforce(outputfn != "", "Output not defined");
  enforce(type == "fasthash", "Algorithm not defined");
  enforce(max_cov != 0, "Maximum coverage not set");
  auto infns = args[1..$];

  assert(max_cov > 0);
  GC.disable();
  foreach (string fn; infns) {
    enforce(outputfn != fn,"Input file can not be same as output file "~fn);
    auto pileup = new PileUp!ReadState();
    auto reader = BamReadBlobStream(fn);
    auto writer = BamWriter(outputfn,reader.header,9);

    // pileup.current = first;

    while(true) { // keep reading while !pileup.empty && !stream.empty
      foreach_invalid_read(reader, (ProcessReadBlob read) {
          writeln("HEY",read);
          auto mod = ModifyProcessReadBlob(read);
          auto blob = mod.toBlob;
          // another hack for now:
          writer.bgzf_writer.write!int(cast(int)(blob.length+2*int.sizeof));
          writer.bgzf_writer.write!int(cast(int)read.raw_ref_id);
          writer.bgzf_writer.write!int(cast(int)read.raw_start_pos);
          writer.bgzf_writer.write(blob);

          // force_write(pileup);
          // current = next(read);
          // pileup.push(current);
        });
      return 1;
    }
  }
  return 0;
}
