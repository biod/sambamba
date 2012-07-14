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
import samheader;
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
    // TODO: check that all headers are valid
    // TODO: check that all files have the same sorting order

    auto ref_id_map = new size_t[size_t][files.length];
    auto program_id_map = new string[string][files.length];
    auto readgroup_id_map = new string[string][files.length];

    auto merged_header = new SamHeader();

    void mergeSequenceDictionaries() {
        static SqLineDictionary mergeTwoDictionaries(SqLineDictionary d1, SqLineDictionary d2)
        {
            auto result = d1;
            int last_index_of_sequence_found_in_d1 = -1; // for checking order inconsistency
           
            foreach (sequence; d2) {
                auto index_in_d1 = result.getSequenceIndex(sequence.name);
                if (index_in_d1 == -1) { // not found
                    result.add(sequence);
                } else {
                    if (index_in_d1 < last_index_of_sequence_found_in_d1) {
                        // can't merge because of issues with sorting order
                        throw new Exception("can't merge SAM headers: reference sequences " ~
                                d1.getSequence(index_in_d1).name ~ " and " ~
                                sequence.name ~ " have inconsistent order");
                    }

                    last_index_of_sequence_found_in_d1 = index_in_d1; // store last index

                    auto expected_sequence_length = d1.getSequence(index_in_d1).length;
                    if (expected_sequence_length != sequence.length) {
                        // those two @SQ lines are highly unlikely to refer to the same
                        // reference sequence if lengths are different
                        throw new Exception("can't merge SAM headers: one of references with " ~
                                "name " ~ sequence.name ~ " has length " ~ 
                                to!string(expected_sequence_length) ~ 
                                " while another one with the same name has length " ~ 
                                to!string(sequence.length));
                    } 
                    // TODO: otherwise, it would be nice to merge individual tags for
                    //       these two @SQ lines
                }
            }

            return result;
        }

        // save new dictionary in merged header
        merged_header.sequences = reduce!mergeTwoDictionaries(map!"a.sequences"(headers));

        // make mapping
        foreach (size_t i, header; headers) {
            foreach (size_t j, sq; header.sequences) {
                auto new_index = merged_header.sequences.getSequenceIndex(sq.name);
                ref_id_map[i][j] = to!size_t(new_index);
            }
        }
    }

    static Tuple!(HeaderLineDictionary!Line, string[string][]) 
    mergeHeaderLines(Line, R)(R records_with_file_ids, size_t file_count)
        if (is(typeof(Line.identifier) == string) &&
            is(ElementType!R == Tuple!(Line, size_t)) &&
            (is(Line == RgLine) || is(Line == PgLine)))
    {
        // Map: record identifier -> record -> list of files
        size_t[][Line][string] id_to_record;

        foreach (record_and_file; records_with_file_ids) {
            auto rec = record_and_file[0];
            auto file_id = record_and_file[1];
            id_to_record[rec.identifier][rec] ~= file_id;
        }

        bool[string] already_used_ids;

        auto record_id_map = new string[string][file_count];
        auto dict = new HeaderLineDictionary!Line;

        // Loop through all identifiers
        foreach (record_id, records_with_same_id; id_to_record) {

            // Several read groups/program records can share the 
            // common identifier, and each one of them can be 
            // presented in several files.
            // 
            // If read groups/program records are equal 
            // (i.e. all fields are equal) then they are treated 
            // as a single read group/program record
            // 
            // Here we iterate over those read groups/program records
            // and files where they were seen, renaming identifiers
            // in order to avoid collisions where necessary.
            foreach (rec, file_ids; records_with_same_id) {
                string new_id = record_id;
                if (record_id !in already_used_ids) {
                    already_used_ids[record_id] = true;
                } else {
                    // if already used ID is encountered,
                    // find unused ID by adding ".N" to the old ID
                    for (int i = 1; ; ++i) {
                        new_id = record_id ~ "." ~ to!string(i);
                        if (new_id !in already_used_ids) {
                            already_used_ids[new_id] = true;
                            break;
                        }
                    }
                }

                // save mapping
                foreach (file_id; file_ids) {
                    record_id_map[file_id][record_id] = new_id;
                }

                // update merged header
                rec.identifier = new_id;
                dict.add(rec);
            }
        }

        return tuple(dict, record_id_map);
    }

    void mergeReadGroups() {
        auto readgroups_with_file_ids = joiner(
            map!((size_t i) {
                    return map!((RgLine line) {
                                   return tuple(line, i);
                                })(headers[i].read_groups.values);
                 })(iota(filenames.length))
        );                             

        auto dict_and_mapping = mergeHeaderLines!RgLine(readgroups_with_file_ids, files.length);
        merged_header.read_groups = dict_and_mapping[0];
        readgroup_id_map = dict_and_mapping[1];
    }

    mergeSequenceDictionaries();
    mergeReadGroups();

    auto alignmentranges_with_file_ids = array(
        map!((size_t i) {
                return tuple(files[i].alignments, i); 
             })(iota(files.length))
    );

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
