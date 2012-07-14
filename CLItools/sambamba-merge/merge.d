/*
   Merging.

   In order for several BAM files to be merged into one, they must firstly be
   sorted in the same order.

   The algorithms used here are based on those from Picard.

   1) Merging headers.
        
        a) Merging sequence dictionaries.
          
            Create map: file -> old ref. ID -> index of sequence name in merged dict

            If sorting order is coordinate, the following invariant holds:
            for any file, reference sequences sorted by old ref. ID appear in 
            the same order in the merged list (if that's impossible, 
            throw an exception because one of files needs to be sorted again
            with a different order of reference sequences)

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
                  -> writeBAM with merged header and reference sequences info

   */

import bamfile;
import bamoutput;

import std.stdio;
import std.algorithm;
import std.conv;
import std.array;
import std.range;
import std.typecons;
import std.stream;

import common.comparators;
import common.nwayunion : nWayUnion;

import utils.samheadermerger;

void printUsage() {
}

int main(string[] args) {

    if (args.length < 3) {
        printUsage();
        return 1;
    }

    auto output_filename = args[1];
    auto filenames = args[2 .. $];
    auto files = array(map!((string fn) { return BamFile(fn); })(filenames));
    auto headers = array(map!"a.header"(files));

    auto merger = new SamHeaderMerger(headers);
    auto merged_header = merger.merged_header;
    auto ref_id_map = merger.ref_id_map;
    auto readgroup_id_map = merger.readgroup_id_map;
    auto program_id_map = merger.program_id_map;

    // tuples of (alignments, file_id)
    auto alignmentranges_with_file_ids = array(
        map!((size_t i) {
                return tuple(files[i].alignments, i); 
             })(iota(files.length))
    );

    // ranges with replaced reference ID, PG and RG tags
    auto modifiedranges = array(
        map!((Tuple!(typeof(BamFile.alignments), size_t) alignments_with_file_id) {
            auto alignments = alignments_with_file_id[0];
            auto file_id = alignments_with_file_id[1];
           
            return map!(
                (Alignment al) {
                    // change reference ID
                    auto old_ref_id = al.ref_id;
                    if (old_ref_id != -1) {
                        auto new_ref_id = to!int(ref_id_map[file_id][old_ref_id]);
                        if (new_ref_id != old_ref_id) {
                            al.ref_id = new_ref_id;
                        }
                    }

                    // change PG tag if it exists
                    auto program = al["PG"];
                    if (!program.is_nothing) {
                        auto pg_str = cast(string)program;
                        auto new_pg = program_id_map[file_id][pg_str];
                        if (new_pg != pg_str) {
                            al["PG"] = new_pg;
                        }
                    }

                    // change reference ID
                    auto read_group = al["RG"];
                    if (!read_group.is_nothing) {
                        auto rg_str = cast(string)read_group;
                        auto new_rg = readgroup_id_map[file_id][rg_str];
                        if (new_rg != rg_str) {
                            al["RG"] = new_rg;
                        }
                    }
                    return al;
                })(alignments);
        })(alignmentranges_with_file_ids)
    );

    // write BAM file

    Stream stream = new BufferedFile(output_filename, FileMode.Out);
    scope(exit) stream.close();

    writeBAM(stream, 
             toSam(merged_header), 
             array(map!((SqLine line) { 
                            ReferenceSequenceInfo ri;
                            ri.name = line.name;
                            ri.length = line.length;
                            return ri;
                        })(merged_header.sequences.values)
             ),
             nWayUnion!compareAlignmentCoordinates(modifiedranges));

    return 0;
}
