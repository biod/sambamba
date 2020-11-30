/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.utils.samheadermerger;

import bio.std.hts.sam.header;
import bio.std.hts.bam.validation.samheader;

import std.array;
import std.range;
import std.algorithm;
import std.conv;
import std.typecons;
import std.exception;

import bio.std.hts.utils.graph;

/// Class encapsulating functionality of merging several SAM headers
/// into one. (In fact, its role is to just group several variables,
/// so it could be replaced by a function returning a struct.)
class SamHeaderMerger {

    enum Strategy {
        simple,
        usingIndex
    }

    Strategy strategy;

    /// Takes array of SAM headers as an input.
    this(SamHeader[] headers, bool validate_headers=false) {
        _headers = headers;
        _len = _headers.length;

        merged_header = new SamHeader();
        ref_id_map = new size_t[size_t][_len];
        ref_id_reverse_map = new size_t[size_t][_len];
        program_id_map = new string[string][_len];
        readgroup_id_map = new string[string][_len];

        if (validate_headers) {
            // TODO: make custom validator for producing better error message
            foreach (size_t i, header; _headers) {
                if (!isValid(header)) {
                    throw new Exception("header #" ~ to!string(i) ~ " is invalid, can't merge");
                }
            }
        }

        auto expected = _headers[0].sorting_order;
        if (expected != SortingOrder.coordinate && expected != SortingOrder.queryname) {
            throw new Exception("file headers indicate that some files are not sorted");
        }
        foreach (header; _headers) {
            if (header.sorting_order != expected) {
                throw new Exception("sorting orders of files don't agree, can't merge");
            }
        }
        merged_header.sorting_order = expected;

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

    /// Map: index of SamHeader in input array of headers -> new refID -> old refID
    size_t[size_t][] ref_id_reverse_map;

private:

    // NOTE FOR DEVELOPER:
    // for more info on what's going on here, read comments in sambamba/sambamba/merge.d

    SamHeader[] _headers;
    size_t _len; // number of headers

    static void addVerticeToDict(ref SqLine[string] dict, ref SqLine line) {
        if (line.name in dict) {
            if (line.length != dict[line.name].length) {
                // those two @SQ lines are highly unlikely to refer to the same
                // reference sequence if lengths are different
                throw new Exception("can't merge SAM headers: one of references with " ~
                        "name " ~ line.name ~ " has length " ~ 
                        to!string(dict[line.name].length) ~ 
                        " while another one with the same name has length " ~ 
                        to!string(line.length));
            }
            // TODO: merge individual tags?
        } else {
            dict[line.name] = line;
        }
    }

    void mergeSequenceDictionaries() {
        // make a directed graph out of reference sequences and do a topological sort
       
        SqLine[string] dict;

        // create a graph
        auto g = new DirectedGraph();
        foreach (header; _headers) {
            auto sequences = header.sequences.values;
            auto prev = sequences.front;
            addVerticeToDict(dict, prev);
            g.addNode(prev.name);
            sequences.popFront();
            while (!sequences.empty) {
                auto cur = sequences.front;
                addVerticeToDict(dict, cur);
                g.addEdge(prev.name, cur.name);
                prev = cur;
                sequences.popFront();
            }
        }

        // get topologically sorted nodes
        try {
            foreach (v; g.topologicalSort()) {
                merged_header.sequences.add(dict[v]);
            }
            strategy = Strategy.simple;
        } catch (Exception e) {
            // failed, try another strategy which requires index files
            foreach (sq_line; sort!((a, b) => a.name < b.name)(dict.values)) {
                merged_header.sequences.add(sq_line);
            }
            strategy = Strategy.usingIndex;
        }

        // make mappings
        foreach (size_t i, header; _headers) {
            foreach (size_t j, SqLine sq; header.sequences) {
                auto new_index = merged_header.sequences.getSequenceIndex(sq.name);
                if (new_index < 0) {
                  import std.stdio;
                  stderr.writeln("merged header sequence dictionary: \n",
                      merged_header.sequences.values);
                  throw new Exception("BUG: " ~ sq.name ~ " is not in merged header dictionary");
                }
                ref_id_map[i][j] = to!size_t(new_index);
                ref_id_reverse_map[i][to!size_t(new_index)] = j;
            }
        }
    }

    // The reason to pass by reference is that when merging program records,
    // this function is called in a loop, and we need to keep some structures between calls.
    //
    // $(D dict) is a collection of Line structs, which will finally be part of the header;
    // $(D record_id_map) is an array of mappings (for each header) where old record identifier
    // is mapped into a new one;
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
                if (record_id in dict) {
                    // if already used ID is encountered,
                    // find unused ID by adding ".N" to the old ID
                    for (int i = 1; ; ++i) {
                        new_id = record_id ~ "." ~ to!string(i);
                        if (new_id !in dict) {
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
        Tuple!(RgLine, size_t)[] readgroups_with_file_ids;
        for (size_t i = 0; i < _len; i++)
            foreach (rg; _headers[i].read_groups.values)
                readgroups_with_file_ids ~= tuple(rg, i);

        auto dict = new RgLineDictionary();

        mergeHeaderLines(readgroups_with_file_ids, _len,
                         dict, readgroup_id_map);

        merged_header.read_groups = dict;
    }

    void mergeProgramRecords() {
        Tuple!(PgLine, size_t)[] programs_with_file_ids;
        for (size_t i = 0; i < _len; i++)
            foreach (pg; _headers[i].programs.values)
                programs_with_file_ids ~= tuple(pg, i);

        auto vertices = partition!"a[0].previous_program !is null"(programs_with_file_ids);
        programs_with_file_ids = programs_with_file_ids[0 .. $ - vertices.length];

        auto dict = new PgLineDictionary();

        while (!vertices.empty) {
            // populates dict and program_id_map
            mergeHeaderLines!PgLine(vertices, _len, dict, program_id_map); 

            // find children of current vertices
            auto old_ids = map!"tuple(a[0].identifier, a[1])"(vertices);
            vertices = partition!((Tuple!(PgLine, size_t) a) {
                                    return !canFind(old_ids, tuple(a[0].previous_program, a[1]));
                                  })(programs_with_file_ids);
            programs_with_file_ids = programs_with_file_ids[0 .. $ - vertices.length];

            // update PP tags in children
            
            foreach (ref pg_with_file_id; vertices) {
                auto pg = pg_with_file_id[0];
                auto file_id = pg_with_file_id[1];
               
                if (pg.previous_program !is null && 
                    pg.previous_program in program_id_map[file_id]) 
                {
                    auto new_id = program_id_map[file_id][pg.previous_program];
                    if (new_id != pg.previous_program) {
                        pg.previous_program = new_id;
                    }
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

unittest {
    import std.stdio;
    import std.algorithm;

    // stderr.writeln("Testing SAM header merging...");
    auto h1 = new SamHeader();
    auto h2 = new SamHeader();
    auto h3 = new SamHeader();
    h1.sorting_order = SortingOrder.coordinate;
    h2.sorting_order = SortingOrder.coordinate;
    h3.sorting_order = SortingOrder.coordinate;
 
    // ---- fill reference sequence dictionaries -------------------------------

    h1.sequences.add(SqLine("A", 100));
    h1.sequences.add(SqLine("B", 200));
    h1.sequences.add(SqLine("C", 300));

    h2.sequences.add(SqLine("D", 100));
    h2.sequences.add(SqLine("B", 200));
    h2.sequences.add(SqLine("E", 300));

    h3.sequences.add(SqLine("A", 100));
    h3.sequences.add(SqLine("E", 300));
    h3.sequences.add(SqLine("C", 300));

    // expected:        A       B       C
    //                      D       E

    // ---- add a few read group records ---------------------------------------

    h1.read_groups.add(RgLine("A", "CN1"));
    h1.read_groups.add(RgLine("C", "CN3"));

    h2.read_groups.add(RgLine("B", "CN2"));
    h2.read_groups.add(RgLine("C", "CN4"));

    h3.read_groups.add(RgLine("B", "CN2"));
    h3.read_groups.add(RgLine("A", "CN4"));

    // ---- add some program records with a lot of name clashes ----------------

    h1.programs.add(PgLine("A", "X"));             //        .> C
    h1.programs.add(PgLine("B", "Y", "", "A"));    //       /
    h1.programs.add(PgLine("C", "Z", "", "B"));    // A -> B -> D
    h1.programs.add(PgLine("D", "T", "", "B"));    //

    h2.programs.add(PgLine("B", "Z"));             // B -> A -> C
    h2.programs.add(PgLine("A", "Y", "", "B"));
    h2.programs.add(PgLine("C", "T", "", "A"));

    h3.programs.add(PgLine("D", "Y"));             // D -> C -> B
    h3.programs.add(PgLine("C", "T", "", "D"));
    h3.programs.add(PgLine("B", "X", "", "C"));

    // expected result:
    //
    //            .> C.1
    //           /
    //  A -> B.1  -> D.1
    //
    //  B -> A.1 -> C.2
    //
    //  D -> C -> B.2

    // ---- add some comments - just for the sake of completeness --------------

    h1.comments ~= "abc";
    h2.comments ~= ["def", "ghi"];
    
    // ------------------ merge these three headers ----------------------------

    {
    auto merger = new SamHeaderMerger([h1, h2, h3]);
    auto h = merger.merged_header;

    assert(h.sorting_order == SortingOrder.coordinate);

    assert(equal(h.sequences.values, 
                 [SqLine("A", 100), SqLine("D", 100), SqLine("B", 200),
                  SqLine("E", 300), SqLine("C", 300)]));

    assert(h.comments == ["abc", "def", "ghi"]);

    assert(equal(sort(array(map!"a.identifier"(h.programs.values))),
                 ["A", "A.1", "B", "B.1", "B.2", "C", "C.1", "C.2", "D", "D.1"]));

    assert(equal(sort(array(map!"a.identifier"(h.read_groups.values))),
                 ["A", "A.1", "B", "C", "C.1"]));
    }

    // sambamba issue 110
    {
    auto h0 = new SamHeader();
    h0.sorting_order = SortingOrder.coordinate;
    h0.sequences.add(SqLine("A", 100));

    auto merger = new SamHeaderMerger([h0]);
    auto h = merger.merged_header;
    assert(equal(h.sequences.values, h0.sequences.values));
    }
}
