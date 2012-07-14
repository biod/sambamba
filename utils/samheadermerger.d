module utils.samheadermerger;

import samheader;

import std.array;
import std.range;
import std.algorithm;
import std.conv;
import std.typecons;

/// Class encapsulating functionality of merging several SAM headers
/// into one. (In fact, its role is to just group several variables,
/// so it could be replaced by a function returning a struct.)
class SamHeaderMerger {

    /// Takes array of SAM headers as an input.
    ///
    /// WARNING: merger might modify the passed array for better performance.
    this(SamHeader[] headers) {
        _headers = headers;
        _len = _headers.length;

        merged_header = new SamHeader();
        ref_id_map = new size_t[size_t][_len];
        program_id_map = new string[string][_len];
        readgroup_id_map = new string[string][_len];

        // TODO: check that all headers are valid
        // TODO: check that all files have the same sorting order

        mergeSequenceDictionaries();
        mergeReadGroups();
        mergeProgramRecords();
        mergeComments();
    }

    /// The main result of merging -- new SamHeader
    SamHeader merged_header;

    /// Map: index of SamHeader in input array of headers -> old refID -> new refID
    size_t[size_t][] ref_id_map;

    /// the same for read group identifiers
    string[string][] readgroup_id_map;

    /// the same for program record identifiers
    string[string][] program_id_map;

private:

    // NOTE FOR DEVELOPER:
    // for more info on what's going on here, read comments in CLItools/sambamba-merge/merge.d

    SamHeader[] _headers;
    size_t _len; // number of headers

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
        merged_header.sequences = reduce!mergeTwoDictionaries(map!"a.sequences"(_headers));

        // make mapping
        foreach (size_t i, header; _headers) {
            foreach (size_t j, sq; header.sequences) {
                auto new_index = merged_header.sequences.getSequenceIndex(sq.name);
                ref_id_map[i][j] = to!size_t(new_index);
            }
        }
    }

    static void mergeHeaderLines(Line, R)(R records_with_file_ids, size_t file_count,
                                          ref HeaderLineDictionary!Line dict,
                                          ref string[string][] record_id_map)
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
    }

    void mergeReadGroups() {
        auto readgroups_with_file_ids = joiner(
            map!((size_t i) {
                    return map!((RgLine line) {
                                   return tuple(line, i);
                                })(_headers[i].read_groups.values);
                 })(iota(_len))
        );                             

        auto dict = new RgLineDictionary();
        mergeHeaderLines(readgroups_with_file_ids, _len,
                         dict, readgroup_id_map);
        merged_header.read_groups = dict;
    }

    void mergeProgramRecords() {
        auto programs_with_file_ids = array(joiner(
            map!((size_t i) {
                    return map!((PgLine line) {
                                   return tuple(line, i);
                                })(_headers[i].programs.values);
                 })(iota(_len))
        ));

        auto vertices = partition!"a[0].previous_program !is null"(programs_with_file_ids);
        programs_with_file_ids = programs_with_file_ids[0 .. $ - vertices.length];

        auto dict = new PgLineDictionary();

        while (!vertices.empty) {
            // populates dict and program_id_map
            mergeHeaderLines!PgLine(vertices, _len, dict, program_id_map); 

            // find children of current vertices
            auto old_ids = map!"a[0].identifier"(vertices);
            vertices = partition!((Tuple!(PgLine, size_t) a) {
                                    return canFind(old_ids, a[0].previous_program);
                                  })(programs_with_file_ids);
            programs_with_file_ids = programs_with_file_ids[0 .. $ - vertices.length];

            // update PP tags in children
            assert(overlap(programs_with_file_ids, vertices));

            foreach (ref pg_with_file_id; vertices) {
                auto pg = pg_with_file_id[0];
                auto file_id = pg_with_file_id[1];
               
                if (pg.previous_program !is null) {
                    pg.previous_program = program_id_map[file_id][pg.previous_program];
                }

                pg_with_file_id = tuple(pg, file_id);
            }
        }

        merged_header.programs = dict;
    }

    void mergeComments() {
        merged_header.comments = join(map!"a.comments"(_headers));
    }
}
