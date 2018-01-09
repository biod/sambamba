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

enum RState { unknown, keep, drop }

// ReadState keeps track of the state of a processed Read. This state
// is maintained on the ringbuffer. We may change this design later.

struct ReadState {
  ProcessReadBlob read;
  RState state;

  this(ProcessReadBlob _r) {
    read = _r;
    state = RState.unknown;
  }

  // @disable this(this); // disable copy semantics;

  @property ref ProcessReadBlob get() {
    return read;
  }

  @property void set_keep() {
    state = RState.keep;
  }
  @property void set_drop() {
    state = RState.drop;
  }
  @property bool is_dropped() {
    return state == RState.drop;
  }
}


/**
   Implementation of Hash1. While reads stream in they get piled up in
   a ringbuffer. The ringbuffer gets filled ahead until the read that is
   no longer overlapping, creating a window:

                                                          r-------------
                                               ----y---------------
                                        -----------y--------
                                 x=================y
                            -----x-------------
                   --------------x----
              l------------------x-----
          leftmost             start_pos       end_pos  rightmost

   once the reads have been counted it moves to the next read, reads ahead
   and destroys the leftmost reads that have gone out of the window.
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
  enforce(type != "", "Algorithm not defined");
  enforce(max_cov != 0, "Maximum coverage not set");
  auto infns = args[1..$];

  assert(max_cov > 0);
  foreach (string fn; infns) {
    enforce(outputfn != fn,"Input file can not be same as output file "~fn);
    auto pileup = new PileUp!ReadState();
    auto stream = BamReadBlobStream(fn);

    auto output = BamWriter(outputfn,stream.header,9);
    auto current = ProcessReadBlob(stream.read);
    asserte(!current.isNull);
    auto current_idx = pileup.push(ReadState(current));
    assert(current_idx == 0);
    // writeln("First pos is ",current.ref_id,":",current.start_pos);
    auto rightmost = current;
    auto rightmost_idx = current_idx;
    auto leftmost = current;
    auto leftmost_idx = current_idx;

    auto reap = () {
      assert(!leftmost.isNull);
      auto readinfo = pileup.read_at_idx(leftmost_idx);
      if (readinfo.is_dropped) {
        writeln("Dropped pos ",leftmost.refid,":",leftmost.start_pos," ",readinfo.state);
      } else {
        if (leftmost.is_mapped)
          writeln("Writing pos ",leftmost.refid,":",leftmost.start_pos," ",readinfo.state);
        else
          writeln("Writing unmapped ",readinfo.state);
      }
      if (!remove) {
        auto mod = ModifyProcessReadBlob(leftmost);
        if (readinfo.is_dropped)
          mod.set_qc_fail;
        auto blob = mod.toBlob;
        // another hack for now:
        output.bgzf_writer.write!int(cast(int)(blob.length+2*int.sizeof));
        output.bgzf_writer.write!int(cast(int)leftmost.raw_ref_id);
        output.bgzf_writer.write!int(cast(int)leftmost.raw_start_pos);
        output.bgzf_writer.write(blob);
      }
      leftmost_idx = pileup.popFront();
      if (!pileup.empty)
        leftmost = pileup.front.get;
    };

    while (true) { // loop through pileup
      assert(!current.isNull);
      while (current.is_unmapped2) {
        // we hit an unmapped set, need to purge (this won't work on threads)
        // writeln("Skip unmapped read");
        reap();
        current = ProcessReadBlob(stream.read);
        if (current.isNull)
          break;
        current_idx = pileup.push(ReadState(current));
        rightmost = current;
        rightmost_idx = current_idx;
      }
      assert(current.is_mapped2);
      // writeln("Current pos is ",current.ref_id,":",current.start_pos);
      // Fill ring buffer ahead until the window is full (current and rightmost)
      // rightmost is null at the end of the genome
      while (!rightmost.isNull && rightmost.is_mapped2 && current.ref_id == rightmost.ref_id && rightmost.start_pos < current.end_pos+1) {
        rightmost = ProcessReadBlob(stream.read);
        if (rightmost.isNull)
          break;
        rightmost_idx = pileup.push(ReadState(rightmost));
      }

      if (false) {
      // Now we have a pileup and we can check this read (output)
      assert(leftmost.is_mapped2);
      writeln("     start  at ",leftmost.ref_id," ",leftmost.start_pos,":",leftmost.end_pos);
      writeln("---> pileup at ",current.ref_id," ",current.start_pos,":",current.end_pos);
      if (!rightmost.isNull && rightmost.is_mapped2)
        writeln("     ending at ",rightmost.ref_id," ",rightmost.start_pos,":",rightmost.end_pos);
      else
        writeln("     reached end ",pileup.ring.length());
      }

      writeln("Current: ",current.show_flags);
      if (!current.is_qc_fail) {
        // Compute depth (leftmost, current, rightmost)
        auto depth = 0;
        auto ldepth = 0;
        auto rdepth = 0;
        for (RingBufferIndex idx = leftmost_idx; idx < rightmost_idx; idx++) {
          auto check = pileup.read_at_idx(idx).get;
          writeln("Check: ",check.show_flags);
          if (check.is_mapped && !check.is_qc_fail) {
            assert(current.is_mapped2);
            assert(check.is_mapped2);
            assert(current.ref_id == check.ref_id);
            // all time is consumed in this section
            // if (reads_overlap(current,check)) { // 8s
            if (read_overlaps(current.start_loc,check)) // 5s
              ldepth++;
            if (read_overlaps(current.end_loc,check)) // 5s
              rdepth++;
            //  depth++;
            // }
          }
        }
        auto this_cov = max(ldepth,rdepth);
        if (this_cov > max_cov) {
          auto hash = SuperFastHash(current.read_name);
          double sample_drop_rate = cast(double)(1 - (this_cov - max_cov)) / this_cov;
          double rand = cast(double)(hash & 0xffffff)/0x1000000;
          write("readinfo ");
          auto readinfo = pileup.read_at_idx(current_idx);
          if (rand < -sample_drop_rate) {
            readinfo.set_drop;
            pileup.update_read_at_index(current_idx,readinfo); // this is ugly
            assert(readinfo.is_dropped);
            write("readinfo2 ");
            auto readinfo2 = pileup.read_at_idx(current_idx);
            writeln(&readinfo, " ", &readinfo2);
            // assert(readinfo is readinfo2);
            assert(readinfo.is_dropped);
            assert(readinfo2.is_dropped);
          }
          else {
            readinfo.set_keep;
            pileup.update_read_at_index(current_idx,readinfo);
          }
          write("Pos ",current.refid,":",current.start_pos);
          writeln(" #",hash," depth ",this_cov," max ",max_cov," sample drop rate ",sample_drop_rate," rand ",rand," ",readinfo.state);
        }
        if (false)
          writeln("**** ",current.read_name," Depth l",ldepth," r",rdepth," t",depth," mapq ",current.mapping_quality()," tlen ", current.tlen," seqlen ",current.sequence_length, " maplen ",current.consumed_reference_bases, " ", current.sequence, "cigar", current.cigar);
      }

      // Stop at end of data
      if (rightmost.isNull && pileup.idx_at_end(current_idx))
        break;

      // Move to next (current)
      current_idx = pileup.get_next_idx(current_idx);
      auto prev = current;
      current = pileup.read_at_idx(current_idx).get;
      if (current.is_mapped2 && prev.is_mapped2 && current.ref_id == prev.ref_id)
        enforce(current.start_pos >= prev.start_pos, "BAM file is not sorted");
      assert(!current.isNull);

      // Reaper: write and remove leading reads (leftmost and current)
      while (!pileup.empty && (leftmost.is_unmapped2 || (leftmost.is_mapped2 && current.is_mapped2 && (leftmost.ref_id != current.ref_id || leftmost.end_pos < current.start_pos)))) {
        // write read
        reap();
      }
    }
    while (!pileup.empty)
      reap();
  }
  return 0;
}

// TODO:
//
//   1. &find template alignment length (end_pos)
//   2. &check depth at &start and &end (should match pileup)
//   3. &quality filter
//   4. &check for valid RNAME in case of CIGAR
//   5. &Write header (bgzf magic), bgzf blocks
//     a. &check ringbuffer implementation
//     b. &create test comparing unpacked versions
//     c. refactor a bit and check for unmapped reads - straighten out flag use
//     d. run memory checker
//   6. Go multi-core
//   7. Introduce option for (development) validation (less checking by default) and
//      introduce assert_throws
//   8. markdup filter
//   9. improve for pairs
