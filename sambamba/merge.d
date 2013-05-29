/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.merge;
/*
   Merging.

   In order for several BAM files to be merged into one, they must firstly be
   sorted in the same order.

   1) Merging headers.
        
        a) Merging sequence dictionaries.
          
            Create map: file -> old ref. ID -> index of sequence name in merged dict

            If sorting order is coordinate, the following invariant holds:
            for any file, reference sequences sorted by old ref. ID appear in 
            the same order in the merged list (if that's impossible, 
            throw an exception because one of files needs to be sorted again
            with a different order of reference sequences)

            For that, a graph is created and topological sort is performed.

        b) Merging program records.
            
            Create map: file -> old program record name -> new program record name

            Some ids can be changed to avoid collisions.

            Invariant: partial order, implied by PP tag, is maintained.
            
            Program records form disjoint set of trees if we consider this relation,
            therefore to maintain partial order we walk these trees with BFS,
            and on each step update changed PP tags.

        c) Merging read group dictionaries.

            The simplest one because there're no restrictions on order.
            Just detect collisions and rename read groups appropriately in such a case.

   2) Merging alignments.

        Use maps built during merging headers:
            file -> old reference id -> new reference id,
            file -> old program record name -> new program record name,
            file -> old read group name -> new read group name

        filenames -> ranges of alignments for these filenames
                  -> ranges of alignments modified accordingly to the maps
                  -> nWayUnion with a comparator corresponding to the common sorting order
                  -> write BAM with merged header and reference sequences info

   */

import bio.bam.reader;
import bio.bam.writer;
import bio.bam.utils.samheadermerger;
import bio.bam.read;

import std.stdio;
import std.algorithm;
import std.conv;
import std.functional;
import std.array;
import std.file;
import std.range;
import std.typecons;
import std.traits;
import std.numeric;
import std.parallelism;
import std.stream;
import std.getopt;

import core.atomic;

import sambamba.utils.common.progressbar;

void printUsage() {
    stderr.writeln("Usage: sambamba-merge [options] <output.bam> <input1.bam> <input2.bam> [...]");
    stderr.writeln();
    stderr.writeln("Options: -t, --nthreads=NTHREADS");
    stderr.writeln("               number of threads to use for compression/decompression");
    stderr.writeln("         -l, --compression-level=COMPRESSION_LEVEL");
    stderr.writeln("               level of compression for merged BAM file, number from 0 to 9");
    stderr.writeln("         -H, --header");
    stderr.writeln("               output merged header to stdout in SAM format, other options are ignored; mainly for debug purposes");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("               show progress bar in STDERR");
}

// these variables can be implicitly used in tasks created in writeBAM
shared(SamHeaderMerger) merger;
shared(SamHeader) merged_header;
shared(size_t[size_t][]) ref_id_map;
shared(string[string][]) program_id_map;
shared(string[string][]) readgroup_id_map;

__gshared static TaskPool task_pool;

BamRead changeAlignment(Tuple!(BamRead, size_t) al_with_file_id) {
    auto al = al_with_file_id[0];
    auto file_id = al_with_file_id[1];
    // change reference ID
    auto old_ref_id = al.ref_id;

    assert(file_id < ref_id_map.length);

    if (old_ref_id != -1 && old_ref_id in ref_id_map[file_id]) {
        auto new_ref_id = to!int(ref_id_map[file_id][old_ref_id]);
        if (new_ref_id != old_ref_id) {
            al.ref_id = new_ref_id;
        }
    } 

    // change PG tag if it exists
    auto program = al["PG"];
    if (!program.is_nothing) {
        auto pg_str = cast(string)program;
        if (pg_str in program_id_map[file_id]) {
            auto new_pg = program_id_map[file_id][pg_str];
            if (new_pg != pg_str) {
                al["PG"] = cast()new_pg;
            }
        }
    }

    // change read group tag
    auto read_group = al["RG"];
    if (!read_group.is_nothing) {
        auto rg_str = cast(string)read_group;
        if (rg_str in readgroup_id_map[file_id]) {
            auto new_rg = readgroup_id_map[file_id][rg_str];
            if (new_rg != rg_str) {
                al["RG"] = cast()new_rg;
            }
        }
    }
    return al;
}

auto modifyAlignmentRange(T)(T alignments_with_file_id) {
    version(serial) {
        return map!changeAlignment(zip(alignments_with_file_id[0], 
                                       repeat(alignments_with_file_id[1])));
    } else {
        return task_pool.map!changeAlignment(zip(alignments_with_file_id[0],
                                                 repeat(alignments_with_file_id[1])),
                                            8192);
    }
}

version(standalone) {
    int main(string[] args) {
        return merge_main(args);
    }
}

int merge_main(string[] args) {

    int compression_level = -1;
    int number_of_threads = totalCPUs;
    bool validate_headers = false;
    bool header_only = false;
    bool show_progress = false;

    if (args.length < 4) {
        printUsage();
        return 1;
    }

    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "nthreads|t",            &number_of_threads,
               "compression-level|l",   &compression_level,
               "validate-headers|v",    &validate_headers,
               "header|H",              &header_only,
               "show-progress|p",       &show_progress);

        task_pool = new TaskPool(number_of_threads);
        scope(exit) task_pool.finish();

        auto output_filename = args[1];
        auto filenames = args[2 .. $];
        BamReader[] files;
        files.length = filenames.length;
        foreach (i; 0 .. files.length) {
            files[i] = new BamReader(filenames[i], task_pool);
            files[i].setBufferSize(50_000_000 / files.length); //TODO
        }
        auto headers = array(map!"a.header"(files));

        merger = cast(shared) new SamHeaderMerger(headers, validate_headers);
        merged_header = merger.merged_header;
        ref_id_map = merger.ref_id_map;
        readgroup_id_map = merger.readgroup_id_map;
        program_id_map = merger.program_id_map;

        if (header_only) {
            write((cast()merged_header).text);
            return 0;
        }

        // alignmentranges must hold tuples of some range of alignment and file id
        void mergeAlignments(R)(R alignmentranges) {
            // ranges with replaced reference ID, PG and RG tags
            auto modifiedranges = array(map!modifyAlignmentRange(alignmentranges));

            // write BAM file
            Stream stream = new BufferedFile(output_filename, FileMode.OutNew, 50_000_000); // TODO
            scope(failure) stream.close();

            auto reference_sequences = new ReferenceSequenceInfo[(cast()merged_header).sequences.length];
            size_t i;
            foreach (line; (cast()merged_header).sequences.values) {
                reference_sequences[i] = ReferenceSequenceInfo(line.name, line.length);
                ++i;
            } 

            auto writer = new BamWriter(stream, compression_level, task_pool);
            scope(exit) writer.finish();
            writer.writeSamHeader(cast()merged_header);
            writer.writeReferenceSequenceInfo(reference_sequences);

            switch (merged_header.sorting_order) {
                case SortingOrder.queryname:
                    foreach (read; nWayUnion!compareReadNames(modifiedranges))
                        writer.writeRecord(read);
                    break;
                case SortingOrder.coordinate:
                    foreach (read; nWayUnion!compareCoordinates(modifiedranges))
                        writer.writeRecord(read);
                    break;
                default: assert(0);
            }
        } // mergeAlignments

        if (show_progress) {
            // tuples of (alignments, file_id)
            shared(float[]) merging_progress;
            merging_progress.length = files.length;

            auto bar = new shared(ProgressBar)();

            alias ReturnType!(BamReader.readsWithProgress!withoutOffsets) AlignmentRangePB;
            auto alignmentranges_with_file_ids = new Tuple!(AlignmentRangePB, size_t)[files.length];

            auto weights = cast(shared)array(map!(pipe!(getSize, to!float))(filenames));
            normalize(cast()weights);

            foreach (i; 0 .. files.length) {
                alignmentranges_with_file_ids[i] = tuple(
                        files[i].readsWithProgress(
                            (size_t j) { 
                                return (lazy float p) {
                                    atomicStore(merging_progress[j], p);
                                    synchronized (bar) {
                                        bar.update(dotProduct(merging_progress, weights));
                                    }
                                };
                            }(i)),
                        i
                );
            }
            mergeAlignments(alignmentranges_with_file_ids);

            bar.finish();
        } else {
            auto alignmentranges_with_file_ids = array(
                zip(map!"a.reads"(files), iota(files.length))
            );
            mergeAlignments(alignmentranges_with_file_ids);
        }
    
    } catch (Throwable e) {
        stderr.writeln("sambamba-merge: ", e.msg);
        return 1;
    }
    return 0;
}
