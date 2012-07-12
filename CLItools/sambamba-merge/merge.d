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
import samheader;

import std.stdio;
import std.algorithm;
import std.array;
import std.range;
import std.typecons;

void printUsage() {
}

int main(string[] args) {

    if (args.length == 1) {
        printUsage();
        return 1;
    }

    auto filenames = args[1 .. $];
    auto files = array(map!((string fn) { return BamFile(fn); })(filenames));
    auto headers = array(map!"a.header"(files));

    auto ref_id_map = new int[int][files.length];
    auto program_id_map = new string[string][files.length];
    auto readgroup_id_map = new string[string][files.length];

    auto merged_header = new SamHeader();

    void mergeReadGroups() {
        auto readgroups_with_file_ids = joiner(
            map!((size_t i) {
                    return map!((RgLine line) {
                                   return tuple(line, i);
                                })(headers[i].read_groups.values);
                 })(iota(filenames.length)));                             
       
        size_t[][RgLine][string] id_to_record;

        foreach (rg_and_file; readgroups_with_file_ids) {
            auto rg = rg_and_file[0];
            auto file_id = rg_and_file[1];
            id_to_record[rg.identifier][rg] ~= file_id;
        }

        bool[string] already_used_ids;
        foreach (rg_id, records_with_same_id; id_to_record) {
            foreach (rg, file_ids; records_with_same_id) {
                string new_id = rg_id;
                if (rg_id !in already_used_ids) {
                    already_used_ids[rg_id] = true;
                } else {
                    for (int i = 1; ; ++i) {
                        new_id = rg_id ~ "." ~ to!string(i);
                        if (new_id !in already_used_ids) {
                            already_used_ids[new_id] = true;
                            break;
                        }
                    }
                }

                foreach (file_id; file_ids) {
                    readgroup_id_map[file_id][rg_id] = new_id;
                }
            }
        }
    }

    mergeReadGroups();
    return 0;
}
