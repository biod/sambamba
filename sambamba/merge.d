/*
    This file is part of Sambamba.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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
import core.memory;

import sambamba.utils.common.progressbar;
import sambamba.utils.common.overwrite;
import sambamba.utils.common.filtering;
import sambamba.utils.common.readstorage;
import bio.core.utils.outbuffer;
import bio.core.utils.roundbuf;
import bio.bam.utils.value;

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
    stderr.writeln("         -F, --filter=FILTER");
    stderr.writeln("               keep only reads that satisfy FILTER");
}

// these variables can be implicitly used in tasks created in writeBAM
shared(SamHeaderMerger) merger;
shared(SamHeader) merged_header;
shared(size_t[size_t][]) ref_id_map;
shared(size_t[size_t][]) ref_id_reverse_map;
shared(string[string][]) program_id_map;
shared(string[string][]) readgroup_id_map;

__gshared static TaskPool task_pool;
__gshared static Filter read_filter;

auto modifyAlignmentRange(T)(T alignments_with_file_id) {
    return AlignmentRangeModifier!T(alignments_with_file_id);
}

ubyte[] modifier(BamRead[] reads, OutBuffer output_buffer, size_t _file_id,
		 ubyte[] tmp)
{
    // it's practically impossible that one will merge > 10000 files,
    // so add 10 extra bytes for each read (5 for RG and 5 for PG)
    size_t bytes_required = 0;
    foreach (read; reads) {
	bytes_required += read.size_in_bytes + 10;
    }

    output_buffer.capacity = bytes_required;

    foreach (read; reads) {
	auto data = read.raw_data;

	auto cigar_offset = 8 * int.sizeof + read.name.length + 1;
	auto tags_offset = cigar_offset + uint.sizeof * read.cigar.length +
	    (3 * read.sequence_length + 1) / 2;

	auto start_offset = output_buffer.data.length + int.sizeof;
	auto p = cast(int*)(output_buffer.data.ptr + start_offset - int.sizeof);
	output_buffer.putUnsafe(int.init);
	size_t chunk_size = read.size_in_bytes - int.sizeof;

	// change ref. id
	auto ref_id_ptr = cast(int*)(data.ptr);
	auto old_ref_id = *ref_id_ptr;
	if (old_ref_id != -1 && old_ref_id in ref_id_map[_file_id]) {
	    auto new_ref_id = to!int(ref_id_map[_file_id][old_ref_id]);
	    if (new_ref_id != old_ref_id)
		*ref_id_ptr = new_ref_id;
	}

	data[cigar_offset - 1] = 0;
	output_buffer.putUnsafe(data[0 .. tags_offset]);

	auto bytes_written = tags_offset;

	// cool, now change PG and RG tags where needed
	foreach (tag, value; read) {
	    Value val = value;
	    if (tag == "PG") {
		auto pg_str = *cast(string*)(&val);
		if (pg_str in program_id_map[_file_id]) {
		    auto new_pg = program_id_map[_file_id][pg_str];
		    if (new_pg != pg_str) {
			val = Value(new_pg);
			auto delta = new_pg.length - pg_str.length;
			chunk_size += cast(int)delta;
		    }
		}
	    } else if (tag == "RG") {
		auto rg_str = *cast(string*)(&val);
		if (rg_str in readgroup_id_map[_file_id]) {
		    auto new_rg = readgroup_id_map[_file_id][rg_str];
		    if (new_rg != rg_str) {
			val = Value(new_rg);
			auto delta = new_rg.length - rg_str.length;
			chunk_size += cast(int)delta;
		    }
		}
	    }
	    output_buffer.putUnsafe(cast(ubyte[])tag);
	    emplaceValue(tmp.ptr, val);
	    auto val_size = sizeInBytes(val);
	    output_buffer.putUnsafe(tmp[0 .. val_size]);
	    bytes_written += 2 + val_size;
	}
	import std.exception;
	enforce(bytes_written == chunk_size);

	*p = cast(int)chunk_size;
    }
    return output_buffer.data;
}

struct AlignmentRangeModifier(T) {
    private {
	InputRange!BamRead _reads;
	size_t _file_id;

	alias TaskWithData!(modifier, size_t, ubyte[]) ModifyTask;
	RoundBuf!ModifyTask _tasks;

	enum _tmp_size = 128 * 1024;
	size_t _task_index;
	ubyte[] _tmp;
	ubyte[] _tmp_space(size_t _task_index) {
	    auto n_tasks = _tmp.length / _tmp_size;
	    auto k =_task_index % n_tasks;
	    return _tmp[k * _tmp_size .. $][0 .. _tmp_size];
	}

	ubyte[] _curr_data;
	size_t _bytes_read;
	size_t _curr_data_len;
	bool _empty;
	BamRead _front;
    }

    bool empty() @property const {
	return _empty;
    }

    this(T alignments_with_file_id) {
	_reads = alignments_with_file_id[0].filtered(read_filter).inputRangeObject;
	_file_id = alignments_with_file_id[1];

	auto n_tasks = max(task_pool.size, 2) * 4;
	_tasks = RoundBuf!ModifyTask(n_tasks);
	_tmp = new ubyte[_tmp_size * n_tasks];

	foreach (i; 0 .. n_tasks) {
	    if (_reads.empty)
		break;
	    auto t = new ModifyTask();
	    t.input_buffer.fill(&_reads);
	    t.run(task_pool, _file_id, _tmp_space(_task_index));
	    _tasks.put(t);
	    ++_task_index;
	}

	popFront();
    }

    BamRead front() {
	return _front;
    }

    void popFront() {
	if (_bytes_read == _curr_data_len) {
	    if (_tasks.empty) {
		_empty = true;
		return;
	    }

	    auto t = _tasks.front;
	    auto data = t.conversion_task.yieldForce();
	    _curr_data.length = max(data.length, _curr_data.length);
	    _curr_data_len = data.length;
	    _curr_data[0 .. data.length] = data[];
	    _bytes_read = 0;
	    _tasks.popFront();
	    if (!_reads.empty) {
		t.input_buffer.clear();
		t.input_buffer.fill(&_reads);
		t.output_buffer.clear();
		t.run(task_pool, _file_id, _tmp_space(_task_index));
		_tasks.put(t);
		++_task_index;
	    }
	}

	int chunk_size = *(cast(int*)(_curr_data.ptr + _bytes_read));
	auto chunk = _curr_data[int.sizeof + _bytes_read .. $][0 .. chunk_size];
	// FIXME big-endian?
	_front = BamRead(chunk);
	_bytes_read += int.sizeof + chunk_size;
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
    string filter_str = null;

    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "nthreads|t",            &number_of_threads,
               "compression-level|l",   &compression_level,
               "validate-headers|v",    &validate_headers,
               "header|H",              &header_only,
               "show-progress|p",       &show_progress,
               "filter|F",              &filter_str);

        if (args.length < 4) {
          printUsage();
          return 1;
        }

        read_filter = createFilterFromQuery(filter_str);

        task_pool = new TaskPool(number_of_threads);
        scope(exit) task_pool.finish();

        auto output_filename = args[1];
        auto filenames = args[2 .. $];
        foreach (filename; filenames)
            protectFromOverwrite(filename, output_filename);

        GC.disable();

        BamReader[] files;
        files.length = filenames.length;
        foreach (i; 0 .. files.length) {
            files[i] = new BamReader(filenames[i], task_pool);
            files[i].setBufferSize(50_000_000 / files.length); //TODO
	    files[i].assumeSequentialProcessing();
        }
        auto headers = array(map!"a.header"(files));

        GC.disable();
        merger = cast(shared) new SamHeaderMerger(headers, validate_headers);
        merged_header = merger.merged_header;
        ref_id_map = merger.ref_id_map;
        ref_id_reverse_map = merger.ref_id_reverse_map;
        readgroup_id_map = merger.readgroup_id_map;
        program_id_map = merger.program_id_map;

        auto strategy = cast()merger.strategy;
        auto n_references = (cast()merged_header).sequences.length;

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

        if (show_progress && strategy == SamHeaderMerger.Strategy.simple) {
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
            if (strategy == SamHeaderMerger.Strategy.simple) {
                auto alignmentranges_with_file_ids = array(
                        zip(map!"a.reads"(files), iota(files.length))
                        );
                mergeAlignments(alignmentranges_with_file_ids);
            } else {
                alias Tuple!(InputRange!BamRead, size_t) R;
                auto alignmentranges_with_file_ids = new R[files.length];
                foreach (k; 0 .. files.length) {
                    size_t[] order; // order of reference ids to fetch
                    order.length = files[k].reference_sequences.length;
                    size_t n_refs; // in this particular file
                    foreach (j; 0 .. n_references) { // refs from all files
                        auto old_ref_id_ptr = j in ref_id_reverse_map[k];
                        if (old_ref_id_ptr !is null) {
                            order[n_refs++] = *old_ref_id_ptr;
                        }
                    }
                    assert(n_refs == order.length);
                    order = order[0 .. n_refs];
                    auto reads = order.map!(
                            ref_id => files[k].reference(cast(int)ref_id)[]
                            ).joiner().chain(files[k].unmappedReads())
                             .map!(r => r.read).inputRangeObject();
                    alignmentranges_with_file_ids[k] = tuple(reads, k);
                }
                mergeAlignments(alignmentranges_with_file_ids);
            }
        }
    
    } catch (Throwable e) {
	version(development) {
	    throw e;
	} else {
	    stderr.writeln("sambamba-merge: ", e.msg);
	    return 1;
	}
    }
    return 0;
}
