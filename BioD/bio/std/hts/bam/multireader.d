/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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
module bio.std.hts.bam.multireader;

import bio.std.hts.sam.header;
import bio.std.hts.bam.reader;
import bio.std.hts.bam.read;
import bio.std.hts.bam.referenceinfo;
import bio.std.hts.utils.samheadermerger;

import std.algorithm;
import std.range;
import std.conv;
import std.parallelism;
import std.array;
import std.numeric : normalize, dotProduct;
import std.exception;
import std.typecons;
import std.file : getSize;

alias size_t FileId;

/// Read from one of multiple BAM files
struct MultiBamRead(R=BamRead) {
    R read;
    alias read this;

    /// from which file it came
    FileId file_id;

    ///
    MultiBamRead dup() @property const {
        return MultiBamRead(read.dup, file_id);
    }
}

// ((MultiBamRead, SamHeaderMerger), (MultiBamRead, SamHeaderMerger)) -> bool
bool compare(T)(auto ref T r1, auto ref T r2) {
    assert(r1[1] == r2[1]);

    SortingOrder sorting_order;
    if (r1[1] is null)
        sorting_order = r1[0].read.reader.header.sorting_order;
    else
        sorting_order = r1[1].merged_header.sorting_order;

    if (sorting_order == SortingOrder.coordinate)
        return compareCoordinatesAndStrand(r1[0], r2[0]);
    else if (sorting_order == SortingOrder.queryname)
        return compareReadNames(r1[0], r2[0]);
    else
        assert(0);
}

// ([BamRead/BamReadBlock], FileId) -> [MultiBamRead]
auto multiBamReads(R)(R reads, FileId index) {
  static if (is(ElementType!R == BamRead))
    return reads.zip(repeat(index)).map!(x => MultiBamRead!BamRead(x[0], x[1]));
  else
    return reads.zip(repeat(index)).map!(x => MultiBamRead!BamRead(x[0].read, x[1]));
}

// (BamReader, SamHeaderMerger, FileId) -> [(MultiBamRead, SamHeaderMerger, FileId)]
auto readRange(BamReader reader, SamHeaderMerger merger, FileId index) {
    return zip(reader.reads.multiBamReads(index), repeat(merger), repeat(index));
}

// (BamReader, SamHeaderMerger, FileId, int, uint, uint) ->
//                                    [(MultiBamRead, SamHeaderMerger, FileId)]
auto readRange(BamReader reader, SamHeaderMerger merger, FileId index,
               int ref_id, uint start, uint end)
{
    int old_ref_id = ref_id;
    if (merger !is null)
        old_ref_id = cast(int)merger.ref_id_reverse_map[index][ref_id];
    auto reads = reader.reference(old_ref_id)[start .. end];
    return zip(reads.multiBamReads(index), repeat(merger), repeat(index));
}

// (BamReader, SamHeaderMerger, FileId, [BamRegion]) ->
//                                    [(MultiBamRead, SamHeaderMerger, FileId)]
auto readRange(BamReader reader, SamHeaderMerger merger, FileId index,
               BamRegion[] regions) {
    if (merger is null) // single reader => no fiddling with ref_id
        return zip(reader.getReadsOverlapping(regions).multiBamReads(index),
                   repeat(merger), repeat(index));

    auto old_regions = merger is null ? regions : regions.dup;
    foreach (j; 0 .. regions.length) {
        auto new_ref_id = regions[j].ref_id;
        if (new_ref_id != -1) {
            auto old_ref_id = cast(uint)merger.ref_id_reverse_map[index][new_ref_id];
            old_regions[j].ref_id = old_ref_id;
        }
    }
    return zip(reader.getReadsOverlapping(old_regions).multiBamReads(index),
               repeat(merger), repeat(index));
}

// ([BamReader], SamHeaderMerger) -> [[(MultiBamRead, SamHeaderMerger, FileId)]]
auto readRanges(BamReader[] readers, SamHeaderMerger merger) {
    return readers.zip(repeat(merger), iota(readers.length))
                  .map!(t => readRange(t[0], t[1], t[2]))();
}

auto readRangeWithProgress
(BamReader reader, SamHeaderMerger merger, FileId index,
 void delegate() f, void delegate(lazy float) u) {
    return zip(reader.readsWithProgress(u, f).multiBamReads(index),
               repeat(merger), repeat(index));
}

auto readRangesWithProgress
(BamReader[] readers, SamHeaderMerger merger,
 void delegate() f, void delegate(lazy float) delegate(size_t) u)
{
    return readers.zip(repeat(merger), iota(readers.length))
                  .map!(t => readRangeWithProgress(t[0], t[1], t[2], f, u(t[2])))();
}

// ([BamReader], SamHeaderMerger, int, uint, uint) ->
//                                    [[(MultiBamRead, SamHeaderMerger, FileId)]]
auto readRanges(BamReader[] readers, SamHeaderMerger merger,
                int ref_id, uint start, uint end)
{
    return readers.zip(repeat(merger), iota(readers.length),
                       repeat(ref_id), repeat(start), repeat(end))
                  .map!(t => readRange(t[0], t[1], t[2], t[3], t[4], t[5]))();
}

// ([BamReader], SamHeaderMerger, [BamRegion]) ->
//        [[(MultiBamRead, SamHeaderMerger, FileId)])
auto readRanges(BamReader[] readers, SamHeaderMerger merger, BamRegion[] regions)
{
    return readers.zip(repeat(merger), iota(readers.length),
                       repeat(regions))
        .map!(t => readRange(t[0], t[1], t[2], t[3]))();
}

// tweaks RG and PG tags, and reference sequence ID
// [[(BamRead, SamHeaderMerger, size_t)]] -> [[MultiBamRead]]
auto adjustTags(R)(R reads_with_aux_info, TaskPool pool, size_t bufsize)
    if (isInputRange!R)
{
  alias R2 = typeof(pool.map!adjustTagsInRange(reads_with_aux_info.front, 1));
  R2[] result;
  foreach (read_range; reads_with_aux_info)
    result ~= pool.map!adjustTagsInRange(read_range, bufsize);
  return result;
}

// (BamRead, SamHeaderMerger, size_t) -> (MultiBamRead, SamHeaderMerger)
auto adjustTagsInRange(R)(R read_with_aux_info) if (!isInputRange!R) {
    auto read = read_with_aux_info[0];
    auto merger = read_with_aux_info[1];
    auto file_id = read_with_aux_info[2];

    if (merger is null) {
        assert(file_id == 0);
        return tuple(read, merger);
    }

    with (merger) {
        assert(file_id < ref_id_map.length);

        auto old_ref_id = read.ref_id;
        if (old_ref_id != -1 && old_ref_id in ref_id_map[file_id]) {
            auto new_ref_id = to!int(ref_id_map[file_id][old_ref_id]);
            if (new_ref_id != old_ref_id)
                read.ref_id = new_ref_id;
        }

        auto program = read["PG"];
        if (!program.is_nothing) {
            auto pg_str = *(cast(string*)(&program));
            if (pg_str in program_id_map[file_id]) {
                auto new_pg = program_id_map[file_id][pg_str];
                if (new_pg != pg_str)
                    read["PG"] = new_pg;
            }
        }

        auto read_group = read["RG"];
        if (!read_group.is_nothing) {
            auto rg_str = *(cast(string*)(&read_group));
            if (rg_str in readgroup_id_map[file_id]) {
                auto new_rg = readgroup_id_map[file_id][rg_str];
                if (new_rg != rg_str)
                    read["RG"] = new_rg;
            }
        }
    }
    return tuple(read, merger);
}

///
class MultiBamReader {

    ///
    this(BamReader[] readers) {
        _readers = readers;

        enforce(_readers.length >= 1, "MultiBamReader requires at least one BAM file");

        if (_readers.length > 1) {
            _merger = new SamHeaderMerger(readers.map!(r => r.header)().array());
            enforce(_merger.strategy == SamHeaderMerger.Strategy.simple, "NYI"); // TODO

            auto n_references = _merger.merged_header.sequences.length;
            _reference_sequences = new ReferenceSequenceInfo[n_references];
            size_t i;
            foreach (line; _merger.merged_header.sequences) {
                _reference_sequences[i] = ReferenceSequenceInfo(line.name, line.length);
                _reference_sequence_dict[line.name] = i++;
            }
        }

        // TODO: maybe try to guess optimal size, based on the number of files?
        setBufferSize(1_048_576);
    }

    ///
    this(string[] filenames) {
        this(filenames.map!(fn => new BamReader(fn))().array());
    }

    ///
    this(string[] filenames, std.parallelism.TaskPool task_pool = taskPool) {
        this(filenames.zip(repeat(task_pool))
                      .map!(fn => new BamReader(fn[0], fn[1]))().array());
    }

    ///
    BamReader[] readers() @property {
        return _readers;
    }

    ///
    SamHeader header() @property {
        return _readers.length > 1 ? _merger.merged_header : _readers[0].header;
    }

    /// Input range of MultiBamRead instances
    auto reads() @property {
        return readRanges(_readers, _merger).adjustTags(task_pool, _adj_bufsz)
                                            .nWayUnion!compare().map!"a[0]"();
    }

    ///
    auto readsWithProgress(void delegate(lazy float p) progressBarFunc,
                           void delegate() finishFunc=null)
    {
        size_t _running = _readers.length;
        void innerFinishFunc() {
            if (--_running == 0 && finishFunc)
                finishFunc();
        }

        auto _progress = new float[_readers.length];
        _progress[] = 0.0;
        auto _weights = _readers.map!(r => r.filename.getSize.to!float).array;
        normalize(_weights);

        auto innerProgressBarFunc(size_t idx) {
            return (lazy float p) {
                _progress[idx] = p;
                progressBarFunc(dotProduct(_progress, _weights));
            };
        }

        return readRangesWithProgress(_readers, _merger,
                                      &innerFinishFunc, &innerProgressBarFunc)
                         .adjustTags(task_pool, _adj_bufsz)
                         .nWayUnion!compare().map!"a[0]"();
    }

    ///
    const(ReferenceSequenceInfo)[] reference_sequences() @property const nothrow {
        if (_readers.length > 1)
            return _reference_sequences;
        else
            return _readers[0].reference_sequences;
    }

    /**
      Check if reference named $(I ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        if (_readers.length > 1)
            return null != (ref_name in _reference_sequence_dict);
        else
            return _readers[0].hasReference(ref_name);
    }

    /**
       Check if all BAM files have indices.
    */
    bool has_index() @property {
        return readers.all!(b => b.has_index);
    }

    /**
      Returns reference sequence with id $(I ref_id).
     */
    MultiBamReference reference(int ref_id) {
        enforce(ref_id >= 0, "Invalid reference index");
        enforce(ref_id < reference_sequences.length, "Invalid reference index");
        return MultiBamReference(_readers, _merger, task_pool, _adj_bufsz,
                                 reference_sequences[ref_id], ref_id);
    }

    /**
      Returns reference sequence named $(I ref_name).
     */
    MultiBamReference opIndex(string ref_name) {
        enforce(hasReference(ref_name),
                "Reference with name " ~ ref_name ~ " does not exist");
        if (_readers.length > 1) {
            auto ref_id = cast(int)_reference_sequence_dict[ref_name];
            return reference(ref_id);
        } else {
            auto ref_id = _readers[0][ref_name].id;
            return reference(ref_id);
        }
    }

    /// Sets buffer size for all readers (default is 1MB)
    void setBufferSize(size_t bytes) {
        foreach (reader; _readers)
            reader.setBufferSize(bytes);
    }

    /**
       Requires coordinate sorting and presence of indices.
    */
    auto getReadsOverlapping(BamRegion[] regions) {
        enforce(header.sorting_order == SortingOrder.coordinate,
                "Not all files are coordinate-sorted");
        enforce(has_index, "Not all files are indexed");

        auto ranges = readRanges(_readers, _merger, regions);
        return ranges.adjustTags(_task_pool, _adj_bufsz)
                     .nWayUnion!compare().map!"a[0]"();
    }

    private {
        BamReader[] _readers;
        SamHeaderMerger _merger;
        ReferenceSequenceInfo[] _reference_sequences;
        size_t[string] _reference_sequence_dict;
        TaskPool _task_pool;
        TaskPool task_pool() @property {
            if (_task_pool is null)
                _task_pool = taskPool;
            return _task_pool;
        }

        size_t _adj_bufsz = 512;
    }
}

///
struct MultiBamReference {
    private {
        BamReader[] _readers;
        SamHeaderMerger _merger;
        int _ref_id;
        ReferenceSequenceInfo _info;
        TaskPool _pool;
        size_t _adj_bufsz;
    }

    this(BamReader[] readers, SamHeaderMerger merger,
         TaskPool task_pool, size_t adj_bufsize,
         ReferenceSequenceInfo info, int ref_id)
    {
        _readers = readers;
        _merger = merger;
        _pool = task_pool;
        _adj_bufsz = adj_bufsize;
        _ref_id = ref_id;
        _info = info;
    }

    ///
    string name() @property const { return _info.name; }

    ///
    int length() @property const { return _info.length; }

    ///
    int id() @property const { return _ref_id; }

    /// Get alignments overlapping [start, end) region.
    /// $(BR)
    /// Coordinates are 0-based.
    auto opSlice(uint start, uint end) {
        enforce(start < end, "start must be less than end");
        auto ranges = readRanges(_readers, _merger, id, start, end);
        return ranges.adjustTags(_pool, _adj_bufsz)
                     .nWayUnion!compare().map!"a[0]"();
    }

    ///
    auto opSlice() {
        return opSlice(0, length);
    }
}
