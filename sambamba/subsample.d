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
  // forwarders
  @property ref ProcessReadBlob get() {
    return read;
  }
  @property GenomePos start_pos() {
    return read.start_pos;
  }
  @property GenomePos end_pos() {
    return read.end_pos;
  }
  @property bool is_mapped() {
    return read.is_mapped;
  }
  @property bool is_unmapped() {
    return read.is_unmapped;
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
     (read+write while is_ignore == unmapped/bad quality/duplicates)
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

void foreach_test_read(ref BamReadBlobStream reader, bool delegate(ProcessReadBlob) test, void delegate(ProcessReadBlob) dg) {
  if (reader.empty) return;
  auto read = ProcessReadBlob(reader.front);
  while(!read.isNull && test(read)) {
    dg(read);
    read = ProcessReadBlob(reader.read);
  }
};

void foreach_test_read(ref BamReadBlobStream reader, bool delegate(ProcessReadBlob) test, bool delegate(ProcessReadBlob) dg) {
  if (reader.empty) return;
  auto read = ProcessReadBlob(reader.front);
  while(!read.isNull && test(read)) {
    if (!dg(read))
      break;
    read = ProcessReadBlob(reader.read);
  }
};

// Same as above but always pop read from reader
void foreach_test_process_read(ref BamReadBlobStream reader, bool delegate(ProcessReadBlob) test, bool delegate(ProcessReadBlob) dg) {
  if (reader.empty) return;
  auto read = ProcessReadBlob(reader.front);
  while(!read.isNull && test(read)) {
    auto res = dg(read);
    read = ProcessReadBlob(reader.read);
    if (!res)
      break;
  }
};

bool do_count(ProcessReadBlob read) {
  return !(read.is_unmapped || read.is_qc_fail || read.is_duplicate);
}

// ---- Move through reads that are ignored.
void foreach_invalid_read(ref BamReadBlobStream reader, void delegate(ProcessReadBlob) dg) {
  foreach_test_read(reader, (read) { return read.is_unmapped || read.is_qc_fail || read.is_duplicate; },dg);
}

// ---- When current is unmapped, move through reads that are ignored.
void foreach_outside_read(ref BamReadBlobStream reader, void delegate(ProcessReadBlob) dg) {
  if (reader.front.is_unmapped) {
    foreach_invalid_read(reader,dg);
  }
}

// ---- Move through reads that are ignored.
void foreach_unmapped_read(ref BamReadBlobStream reader, void delegate(ProcessReadBlob) dg) {
  foreach_test_read(reader, (read) { return read.is_unmapped; },dg);
}

// ---- Move through reads that are ignored.
void foreach_mapped_read(ref BamReadBlobStream reader, void delegate(ProcessReadBlob) dg) {
  foreach_test_read(reader, (read) { return read.is_mapped; },dg);
}

// ---- Read and process each read
void foreach_read(ref BamReadBlobStream reader, bool delegate(ProcessReadBlob) dg) {
  foreach_test_process_read(reader, (read) { return true; },dg);
}


// keep reading until rightmost outside the current read (stack) or until we hit
// a new ref_id
bool keep_reading_while_in_window(PileUp!ReadState pileup, ProcessReadBlob read) {
  auto rightmost = read;
  auto current = pileup.read_current.read;
  asserte(current.is_mapped,"Trying to process an unmapped read");
  if (rightmost.is_unmapped)
    return true; // we have to push on until we have a mapped one
  return current.ref_id == rightmost.ref_id && rightmost.start_pos < current.end_pos+1;
}

struct DepthPos {
  GenomePos pos;
  ulong depth;
}

// Only start_pos depth is stored. If necessary, the end_pos depth is estimated
// from the nearest start positions

class Depth {
  RingBuffer!DepthPos cache; // start_pos is always sorted!
  Nullable!RefId ref_id;

  this(ulong bufsize=DEFAULT_BUFFER_SIZE) {
    cache = RingBuffer!DepthPos(bufsize);
  }

  // Adds a new read, updates depth at read.start and read.end. We are assuming
  // reads coming in are sorted with a valid position
  void add(ProcessReadBlob read) {
    // writeln("Adding ",read.start_pos);
    if (ref_id.isNull || ref_id.get != read.ref_id) {
      // moving into a new window/pileup
      cache.cleanup;
      ref_id = read.ref_id;
    }
    if (cache.empty || cache.back.pos != read.start_pos) {
      // new position to record in the cache
      cache.put(DepthPos(read.start_pos,1));
    }
    else {
      // update existing position
      cache.back.depth++;
    }
    writeln(cache.toString);
  }

  void cleanup_before(GenomePos pos) {
    // pop items of front
    if (cache.empty) return;
    for (DepthPos d = cache.front; d.pos < pos; d = cache.front) {
      cache.popFront;
      if (cache.empty) return;
    }
  }
};

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
  uint count_pileup_full = 0;

  foreach (string fn; infns) {
    enforce(outputfn != fn,"Input file can not be same as output file "~fn);
    auto pileup = new PileUp!ReadState(max_cov * 100);
    auto depth = new Depth(max_cov * 100);
    auto reader = BamReadBlobStream(fn);
    reader.popFront;
    auto writer = BamWriter(outputfn,reader.header,9);

    // The main loop moves pileup.current a step at a time. State is maintained
    // in the reader (current file pos), pileup (leftmost, current, rightmost)
    // and the Depth cache (position depth).
    while(!reader.empty) {
      // Stage1: read ahead for mapped reads and calculate depth
      foreach_read(reader, (ProcessReadBlob read) {
          writeln("Readahead ",read);
          pileup.push(ReadState(read)); // push onto pileup
          if (pileup.is_full)
            return false; // move on for processing
          if (pileup.read_current.is_unmapped)
            return false; // move on to next current read in pileup
          if (do_count(read)) // valid read
            depth.add(read);
          if (keep_reading_while_in_window(pileup,read))
            return true; // keep cycling in loop
          else {
            // moved out of the window
            writeln("Out of window");
            return false;       // and move on for processing
          }
        });

      // Stage2: check pileup buffer is full
      if (pileup.is_full) {
        count_pileup_full++;
        pileup.purge( (ReadState read) {
            write("f");
            writer.push(read.read);
            enforce(false,"Pileup buffer is full");
          });
      }
      // Stage3: mark reads in pileup - there should be no IO here
      if (!pileup.empty) {
        auto read = pileup.read_current;
        if(!pileup.current_is_tail)
          pileup.current_inc;
      }

      // Stage4: write out-of-scope reads and remove from ringbuffer
      pileup.each_left_of_current( (RingBufferIndex idx, ReadState read) {
          read.set_dirty;
          pileup.ring.update_at(idx,read); // because of copy semantics
        });

      pileup.purge_while( (ReadState read) {
          if (read.is_dirty) {
            write("d");
            writer.push(read.read);
            return true;
          }
          return false; // skip rest and do not reset current
        });

      // purge unused positions in Depth ring buffer
      if (!pileup.empty && !pileup.current.isNull && do_count(pileup.read_current.read))
        depth.cleanup_before(pileup.read_current.start_pos); // discard out of window
    }
    // Finally write out remaining reads
    pileup.purge( (ReadState read) {
        write(",");
        writer.push(read.read);
      });
    stderr.writeln(pileup.stats);
  }
  stderr.writeln("Pileup was full ",count_pileup_full," times");
  return 0;
}
