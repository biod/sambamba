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
import std.parallelism;
import std.stream;

import common.comparators;
import common.nwayunion : nWayUnion;

import utils.samheadermerger;

void printUsage() {
    stderr.writeln("Usage: sambamba-merge <output.bam> <input1.bam> [<input.bam>..]");
}

// these variables can be implicitly used in tasks created in writeBAM
shared(SamHeaderMerger) merger;
shared(SamHeader) merged_header;
shared(size_t[size_t][]) ref_id_map;
shared(string[string][]) program_id_map;
shared(string[string][]) readgroup_id_map;

__gshared static TaskPool task_pool;

Alignment changeAlignment(Tuple!(Alignment, size_t) al_with_file_id) {
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

auto modifyAlignmentRange(Tuple!(typeof(BamFile.alignments), size_t) alignments_with_file_id) {
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

    if (args.length < 3) {
        printUsage();
        return 1;
    }

    try {

    task_pool = new TaskPool(totalCPUs);
    scope(exit) task_pool.finish();

    auto output_filename = args[1];
    auto filenames = args[2 .. $];
    BamFile[] files;
    files.length = filenames.length;
    foreach (i; 0 .. files.length) {
        files[i] = BamFile(filenames[i], task_pool);
        files[i].setBufferSize(50_000_000 / files.length); //TODO
    }
    auto headers = array(map!"a.header"(files));

    merger = new shared(SamHeaderMerger)(headers);
    merged_header = merger.merged_header;
    ref_id_map = merger.ref_id_map;
    readgroup_id_map = merger.readgroup_id_map;
    program_id_map = merger.program_id_map;

    // tuples of (alignments, file_id)
    auto alignmentranges_with_file_ids = array(
        zip(map!"a.alignments"(files), iota(files.length))
    );

    // ranges with replaced reference ID, PG and RG tags
    auto modifiedranges = array(map!modifyAlignmentRange(alignmentranges_with_file_ids));

    // write BAM file

    Stream stream = new BufferedFile(output_filename, FileMode.Out, 50_000_000); // TODO
    scope(exit) stream.close();

    auto reference_sequences = new ReferenceSequenceInfo[(cast()merged_header).sequences.length];
    size_t i;
    foreach (line; (cast()merged_header).sequences.values) {
        reference_sequences[i].name = line.name;
        reference_sequences[i].length = line.length;
        ++i;
    } 

    writeBAM(stream, 
             toSam(cast()merged_header), 
             reference_sequences,
             nWayUnion!compareAlignmentCoordinates(modifiedranges),
             -1,
             task_pool);

    } catch (Throwable e) {
        stderr.writeln("sambamba-merge: ", e.msg);
        return 1;
    }
    return 0;
}
