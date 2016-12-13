/*
    This file is part of Sambamba.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.markdup;

import std.stdio;
import std.getopt;
import std.path : buildPath;
import sambamba.utils.common.tmpdir;
import sambamba.utils.common.progressbar;
import sambamba.utils.common.overwrite;
import thirdparty.unstablesort;
import utils.lz4;

import bio.bam.reader, bio.bam.readrange, bio.bam.writer, bio.bam.referenceinfo,
       bio.bam.read, bio.sam.header, bio.bam.abstractreader,
       bio.bam.multireader;
import std.traits, std.typecons, std.range, std.algorithm, std.parallelism, 
       std.exception, std.file, std.typetuple, std.conv, std.array, std.bitmanip,
       std.c.stdlib, std.datetime, std.stream : BufferedFile, FileMode;

/// Read + its index (0-based)
struct IndexedBamRead {
    ulong index;
    BamRead read;
    alias read this;

    IndexedBamRead dup() @property const {
        return IndexedBamRead(index, read.dup);
    }
}

auto withIndices(R)(R reads) {
    return reads.zip(sequence!((a,n)=>n)())
                .map!(t => IndexedBamRead(t[1], t[0]))();
}

struct ReadPair {
    IndexedBamRead read1;
    IndexedBamRead read2;
}

struct ReadPairOrFragment {
    IndexedBamRead read1;
    Nullable!IndexedBamRead read2;

    this(IndexedBamRead r1) {
        read1 = r1;
    }

    this(IndexedBamRead r1, IndexedBamRead r2) {
        read1 = r1;
        read2 = r2;
    }
}

class MallocArray(T) {
    private {
        T* _p;
        size_t _sz;
        size_t _cur;
        float _gf;
    }

    /// call manually!
    void free() {
        std.c.stdlib.free(_p);
        _p = null;
    }

    this(size_t initial_size=16, float grow_factor=1.5) {
        enforce(initial_size >= 16);
        enforce(grow_factor >= 1.2);
        _p = cast(T*)std.c.stdlib.malloc(initial_size * T.sizeof);
        _sz = initial_size;
        _gf = grow_factor;
    }

    T[] data() @property {
        return _p[0 .. _cur];
    }

    size_t capacity() @property const {
        return _sz;
    }

    /// sets length to zero without releasing the memory
    void reset() {
        _cur = 0;
    }

    void put(ref T element) {
        if (_cur < _sz) {
            _p[_cur++] = element;
        } else {
            assert(_cur == _sz);
            _sz = cast(size_t)(_sz * _gf);
            _p = cast(T*)std.c.stdlib.realloc(_p, _sz * T.sizeof);
            _p[_cur++] = element;
        }
    }

    void put(T element) {
        put(element);
    }
}

/// 48 bytes; fast access to read group;
/// also stores precomputed hash.
struct HReadBlock {
    IndexedBamRead read;
    alias read this;
    alias read get;
    uint hash;
    mixin(bitfields!(ushort, "_rg_pos", 16,
                     ushort, "_rg_len", 15,
                     bool, "_is_not_null", 1));
    bool isNull() @property const { return !_is_not_null; }
    string read_group() @property const {
        if (_rg_pos == 0)
            return "";
        return cast(string)read.raw_data[_rg_pos .. _rg_pos + _rg_len];
    }
    void nullify() {
        .destroy(read);
        _is_not_null = false;
    }

    HReadBlock dup() @property {
        auto result = this;
        result.read = read.dup;
        return result;
    }
}
static assert(HReadBlock.sizeof == 40);

template makeHReadBlock(alias charsHashFunc) {
    HReadBlock makeHReadBlock(R)(auto ref R read) {
        HReadBlock r = void;
        r.read = read;
        auto rg = cast(ubyte[])getRG(read);
        if (rg.length > 0) {
            assert(rg.length <= ushort.max / 2);
            assert(rg.ptr - read.raw_data.ptr <= ushort.max);
            r._rg_pos = cast(ushort)(rg.ptr - read.raw_data.ptr);
            r._rg_len = cast(ushort)(rg.length);
        } else {
            r._rg_pos = 0;
        }
        r._is_not_null = true;
        assert(r.read_group == getRG(r.read));
        auto h = charsHashFunc(chain(cast(ubyte[])r.name,
                                     cast(ubyte[])r.read_group));
        r.hash = h & 0xFFFFFFFF;
        return r;
    }
}

private string getRG(R)(auto ref R r) {
    static if (is(R == HReadBlock)) {
        assert(r.read_group == getRG(r.read));
        return r.read_group;
    } else {
        auto v = r["RG"];
        if (v.is_nothing)
            return "";
        return *(cast(string*)(&v));
    }
}

bool compareReadNamesAndReadGroups(R1, R2)(auto ref R1 r1, auto ref R2 r2) {
    auto name_cmp_result = std.algorithm.cmp(r1.name, r2.name);
    if (name_cmp_result < 0)
        return true;
    if (name_cmp_result > 0)
        return false;
    return getRG(r1) < getRG(r2);
}

// LDC doesn't like lambdas :-(
auto _rse2brb(T)(T t) { return IndexedBamRead(t[1], t[0]); }
auto _rse2rs(T)(T t) { return map!_rse2brb(zip(t[0].reads, t[1])); }
auto _rse2r(BamReader[] r, SimpleReader!(ulong, std.stdio.File)[] s) {
    return map!_rse2rs(zip(r, s)).array()
              .nWayUnion!compareReadNamesAndReadGroups();
}

private auto readsFromTempFiles(size_t buf_size, string[] tmp_filenames,
                                TaskPool pool) {
    BamReader[] readers;
    SimpleReader!(ulong, std.stdio.File)[] indices;
    foreach (fn; tmp_filenames) {
        readers ~= new BamReader(fn, pool);
        indices ~= new SimpleReader!(ulong, std.stdio.File)(fn ~ ".idx");
        readers[$ - 1].setBufferSize(buf_size / tmp_filenames.length);
    }

    return _rse2r(readers, indices);
}

struct CollateReadPairRange(R, bool keepFragments, alias charsHashFunc)
    if (isInputRange!R && is(Unqual!(ElementType!R) == IndexedBamRead))
{
    private {
        static auto wrapper(R reads, TaskPool pool) {
            auto r1 = reads.filter!q{ a.ref_id != -1 }
                           .filter!q{ !a.is_unmapped }
                           .filter!q{ !a.is_secondary_alignment }
                           .filter!q{ !a.is_supplementary };

            static if (keepFragments) {
                auto r2 = r1;
            } else {
                auto r2 = r1.filter!q{a.is_paired && !a.mate_is_unmapped}();
            }

            return pool.map!(makeHReadBlock!charsHashFunc)(r2, 1024);
        }

        ReturnType!wrapper _reads;
        TaskPool _task_pool;

        alias IndexedBamRead Read;

        HReadBlock[] _table; // size is always power of two
        size_t _table_mask;

        ulong _min_idx, _max_idx;

        version(profile) {
            size_t _duped_during_compaction;
            StopWatch _compact_sw;
        }

        void _compact() {
            enum max_diff = 2_000UL;
            if (_max_idx - _min_idx < max_diff * 6)
                return;

            version(profile) { _compact_sw.start(); scope(exit) _compact_sw.stop(); }
            
            _min_idx = _max_idx;
            foreach (ref r; _table)
                if (!r.isNull && r.is_slice_backed
                    && _max_idx - r.index > max_diff)
                {
                    version(profile) ++_duped_during_compaction;
                    r = r.dup;
                } else if (!r.isNull && r.is_slice_backed) {
                    _min_idx = min(_min_idx, r.index);
                }
        }

        version(profile) {
            ~this() {
                stderr.writeln("duped during compaction:     ", _duped_during_compaction);
                stderr.writeln("time spent on compaction:    ", _compact_sw.peek().msecs, " ms");
            }
        }

        HReadBlock[] _overflow_list;
        size_t _stored_in_overflow_list; // current number of elements
        size_t _overflow_list_cur_pos; // for finishing current list
        HReadBlock _process_after_dumping;

        enum Source { hashTable, overflowList, tempFiles, none }
        Source _src;

        IBamSamReader _reader;

        static if (keepFragments)
            alias ReadPairOrFragment FrontType;
        else
            alias ReadPair FrontType;
        FrontType _front;

        string _tmp_dir;
        string[] _tmp_filenames;
        BamWriter _tmp_w;
        size_t _tmp_written;

        std.stdio.File _tmp_w_idx;
        MallocArray!ulong _tmp_idx;

        Nullable!Read _tmp_r1;
        ReturnType!readsFromTempFiles _tmp_reads;
    }

    this(R reads, ubyte table_size_log2, size_t overflow_list_size,
         string tmp_dir, TaskPool task_pool = taskPool)
    {
        enforce(overflow_list_size > 0);
        _tmp_dir = tmp_dir;

        _reads = wrapper(reads, task_pool);
        _task_pool = task_pool;
        setSource(Source.hashTable);

        _table = new HReadBlock[1 << table_size_log2];
        _table_mask = (1 << table_size_log2) - 1;
        _overflow_list = new HReadBlock[overflow_list_size];

        _tmp_idx = new MallocArray!ulong(8192);

        popFront();
    }

    bool empty() const { return _src == Source.none; }

    auto front() { return _front; }

    void popFront() {
        final switch (_src) {
            case Source.hashTable:    popFrontHashTable();    return;
            case Source.overflowList: popFrontOverflowList(); return;
            case Source.tempFiles:    popFrontTempFiles();    return;
            case Source.none:         assert(false);
        }
    }

    private:
    void setSource(Source source) {
        _src = source;
    }

    auto next(R)(ref R range) {
        auto r = range.front;
        range.popFront();
        if (_src == Source.tempFiles)
            r.associateWithReader(_reader);
        return r;
    }

    static bool readsArePaired(R1, R2)(auto ref R1 r1, auto ref R2 r2) {
        return r1.name == r2.name && getRG(r1) == getRG(r2);
    }

    size_t computeHash(ref HReadBlock r) {
        return r.hash & _table_mask;
    }

    void copyToOverflowList(ref HReadBlock r) {
        auto old = r.is_slice_backed ? r.dup : r;
        _overflow_list[_stored_in_overflow_list++] = old;
    }

    void updateHashTableEntry(size_t position, ref HReadBlock read) {
        _max_idx = max(_max_idx, read.index);
        _table[position] = read;
    }

    void removeHashTableEntry(size_t position) {
        _table[position].nullify();
    }

    void createTmpWriter() {
        if (_reader is null) {
            _tmp_w = null;
            return;
        }

        _tmp_filenames ~= _tmp_dir ~ "/sorted." ~ 
                          _tmp_filenames.length.to!string() ~ ".bam";
        _tmp_w = new BamWriter(_tmp_filenames[$ - 1], 1, _task_pool);
        _tmp_w_idx = std.stdio.File(_tmp_filenames[$ - 1] ~ ".idx", "w+");
        _tmp_w.disableAutoIndexCreation();
        _tmp_w.writeSamHeader(_reader.header);
        _tmp_w.writeReferenceSequenceInfo(_reader.reference_sequences);
        _tmp_written = 0;
    }

    void closeTmpWriter() {
        if (_tmp_w !is null)
            _tmp_w.finish();
        if (_tmp_filenames.length > 0) {
            _tmp_w_idx.rawWrite(_tmp_idx.data);
            _tmp_w_idx.close();
        }
        _tmp_idx.reset();
    }

    void dumpTmpRecord(R)(auto ref R read) {
        assert(_tmp_w !is null);
        _tmp_w.writeRecord(read);
        _tmp_written++;
        _tmp_idx.put(read.index);
        if (_tmp_idx.data.length == _tmp_idx.capacity) {
            _tmp_w_idx.rawWrite(_tmp_idx.data);
            _tmp_idx.reset();
        }
    }
        
    void popFrontHashTable() {
        while (!_reads.empty) {
            _compact();
            auto read = next(_reads);

            if (_reader is null) {
                _reader = read.reader;
            }

            static if (keepFragments) {
                if (!read.is_paired || read.mate_is_unmapped) {
                    _front = FrontType(read);
                    return;
                }
            }

            auto h = computeHash(read);
            if (_table[h].isNull) {
                updateHashTableEntry(h, read);
            } else if (readsArePaired(_table[h], read)) {
                _front = FrontType(_table[h], read);
                removeHashTableEntry(h);
                return;
            } else if (_stored_in_overflow_list < _overflow_list.length) {
                copyToOverflowList(_table[h]);
                updateHashTableEntry(h, read);
            } else {
                _process_after_dumping = read.dup;
                sort!compareReadNamesAndReadGroups(_overflow_list[]);
                createTmpWriter();
                _overflow_list_cur_pos = 0;
                setSource(Source.overflowList);
                popFrontOverflowList();
                return;
            }
        }

        auto remaining_reads = _table.filter!"!a.isNull"
                                     .array().sort!compareReadNamesAndReadGroups();
        createTmpWriter();
        foreach (r; remaining_reads)
            dumpTmpRecord(r);
        closeTmpWriter();

        if (_stored_in_overflow_list > 0) {
            createTmpWriter();
            auto list = _overflow_list[0 .. _stored_in_overflow_list];
            sort!compareReadNamesAndReadGroups(list);
            foreach (r; list) dumpTmpRecord(r);
            closeTmpWriter();
        }

        _table = null;
        _overflow_list = null;

        // FIXME: constant!
        _tmp_reads = readsFromTempFiles(128_000_000, _tmp_filenames, _task_pool);
        setSource(Source.tempFiles);
        if (!_tmp_reads.empty)
            _tmp_r1 = next(_tmp_reads);
        popFrontTempFiles();
    }

    void popFrontOverflowList() {
        while (_overflow_list_cur_pos < _stored_in_overflow_list - 1) {
            auto r1 = _overflow_list[_overflow_list_cur_pos];
            auto r2 = _overflow_list[_overflow_list_cur_pos + 1];
            if (readsArePaired(r1, r2)) {
                _overflow_list_cur_pos += 2;
                _front = FrontType(r1, r2);
                return;
            } else {
                dumpTmpRecord(_overflow_list[_overflow_list_cur_pos++]);
            }
        }

        if (_overflow_list_cur_pos < _stored_in_overflow_list)
            dumpTmpRecord(_overflow_list[_overflow_list_cur_pos]);
        closeTmpWriter();

        _stored_in_overflow_list = 0;
        if (!_process_after_dumping.isNull) {
            auto h = computeHash(_process_after_dumping);
            copyToOverflowList(_table[h]);
            updateHashTableEntry(h, _process_after_dumping);
            _process_after_dumping.nullify();
        }

        setSource(Source.hashTable);
        popFrontHashTable();
    }

    void popFrontTempFiles() {
        while (!_tmp_reads.empty) {
            auto r2 = next(_tmp_reads);
            if (readsArePaired(_tmp_r1, r2)) {
                _front = FrontType(_tmp_r1, r2);
                if (!_tmp_reads.empty)
                    _tmp_r1 = next(_tmp_reads);
                else
                    _tmp_r1.nullify();
                return;
            } else {
                static if (keepFragments) {
                    _front = FrontType(_tmp_r1);
                    _tmp_r1 = r2;
                    return;
                } else {
                    _tmp_r1 = r2;
                }
            }
        }

        if (!_tmp_r1.isNull) {
            _front = FrontType(_tmp_r1);
            _tmp_r1.nullify();
            return;
        }

        foreach (fn; _tmp_filenames) {
            std.file.remove(fn);
            std.file.remove(fn ~ ".idx");
        }
        setSource(Source.none);
    } 
}

auto simpleHash(R)(R chars) 
if (isInputRange!R && is(ElementType!R == ubyte))
{
    hash_t h = 0;
    foreach (ubyte c; chars) {
        h += c;
        h *= 37;
    }
    return h;
}

auto readPairs(alias hashFunc=simpleHash, R)
(R reads, ubyte table_size_log2, size_t overflow_list_size,
 string tmp_dir, TaskPool task_pool=taskPool) {
    return CollateReadPairRange!(R, false, hashFunc)
            (reads, table_size_log2, overflow_list_size, tmp_dir, task_pool);
}

auto readPairsAndFragments(alias hashFunc=simpleHash, R)
(R reads, ubyte table_size_log2, size_t overflow_list_size,
 string tmp_dir, TaskPool task_pool=taskPool) {
    return CollateReadPairRange!(R, true, hashFunc)
            (reads, table_size_log2, overflow_list_size, tmp_dir, task_pool);
}

/////////////////////////////////////////////////////////////////////////////////
/// no more than 32767 libraries and 16383 reference sequences are supported ////
/////////////////////////////////////////////////////////////////////////////////

// 8 bytes
struct SingleEndBasicInfo {
    mixin(bitfields!(short, "library_id", 16,
                     ushort, "ref_id", 14,
                     ubyte, "reversed", 1,
                     ubyte, "paired", 1));
    int coord;

    bool samePosition(SingleEndBasicInfo other) const {
        return coord == other.coord && ref_id == other.ref_id &&
            reversed == other.reversed && library_id == other.library_id;
    }
}
static assert(SingleEndBasicInfo.sizeof == 8);

// 24 bytes :-(
struct SingleEndInfo {
    SingleEndBasicInfo basic_info;
    alias basic_info this;
    ulong idx;
    ushort score;
}

// 32 bytes
struct PairedEndsInfo {
    mixin(bitfields!(short, "library_id", 16,
                     ushort, "ref_id1", 14,
                     ubyte, "reversed1", 1,
                     ubyte, "reversed2", 1));
    int coord1;
    int coord2;
    ushort ref_id2;

    ushort score; // sum of base qualities that are >= 15 
    ulong idx1, idx2;

    SingleEndBasicInfo read1_basic_info() @property {
        typeof(return) result = void;
        // HACK! HACK! HACK! use the fact that the structures are almost
        // the same except the one bit meaning 'paired' instead of 'reversed2'.
        auto p = cast(ubyte*)(&result);
        p[0 .. 8] = (cast(ubyte*)(&this))[0 .. 8];
        result.paired = true;
        return result;
    }
}

static assert(PairedEndsInfo.sizeof == 32);

bool singleEndInfoComparator(S1, S2)(auto ref S1 s1, auto ref S2 s2) {
    if (s1.library_id < s2.library_id) return true;
    if (s1.library_id > s2.library_id) return false;
    if (s1.ref_id < s2.ref_id) return true;
    if (s1.ref_id > s2.ref_id) return false;
    if (s1.coord < s2.coord) return true;
    if (s1.coord > s2.coord) return false;
    if (s1.reversed < s2.reversed) return true;
    return false;
}

bool pairedEndsInfoComparator(P1, P2)(auto ref P1 p1, auto ref P2 p2) {
    if (p1.library_id < p2.library_id) return true;
    if (p1.library_id > p2.library_id) return false;
    if (p1.ref_id1 < p2.ref_id1) return true;
    if (p1.ref_id1 > p2.ref_id1) return false;
    if (p1.coord1 < p2.coord1) return true;
    if (p1.coord1 > p2.coord1) return false;
    if (p1.reversed1 < p2.reversed1) return true;
    if (p1.reversed1 > p2.reversed1) return false;
    if (p1.reversed2 < p2.reversed2) return true;
    if (p1.reversed2 > p2.reversed2) return false;
    if (p1.ref_id2 < p2.ref_id2) return true;
    if (p1.ref_id2 > p2.ref_id2) return false;
    if (p1.coord2 < p2.coord2) return true;
    return false;
}

struct MarkDuplicatesConfig {
    ubyte hash_table_size_log2 = 18;
    size_t overflow_list_size = 200_000;
    string tmpdir;
    size_t bufsize = 500_000_000; // for sorting storage

    // called on each group of PE duplicates
    void delegate(PairedEndsInfo[]) pe_callback = null;

    // called on each group of SE duplicates;
    // the second argument tells if there is a paired read at that position
    void delegate(SingleEndInfo[], bool) se_callback = null;
}

private /* algorithm */ {

class ReadGroupIndex {
    private {
        int[string] _rg_ids;
        int[] _rg_id_to_lb_id;
    }

    this(SamHeader header) {
        int[string] libraries;
        _rg_id_to_lb_id.length = header.read_groups.length;
        int i = 0;
        int n_libs = 0;
        foreach (rg; header.read_groups) {
            _rg_ids[rg.identifier] = i;
            if (rg.library in libraries) {
                _rg_id_to_lb_id[i] = libraries[rg.library];
            } else {
                _rg_id_to_lb_id[i] = n_libs;
                libraries[rg.library] = n_libs;
                ++n_libs;
            }
            ++i;
        }
        enforce(n_libs <= 32767, "More than 32767 libraries are unsupported");
    }

    /// -1 if read group with such name is not found in the header
    int getId(string name) const {
        auto p = name in _rg_ids;
        if (p is null)
            return -1;
        return *p;
    }

    int getLibraryId(int read_group_id) const {
        if (read_group_id == -1)
            return -1;
        return _rg_id_to_lb_id[read_group_id];
    }
}

int computeFivePrimeCoord(R)(auto ref R read) {
    if (!read.is_reverse_strand) {
        auto ops = read.cigar.until!q{ !a.is_clipping };
        return read.position - reduce!q{ a + b }(0, ops.map!q{ a.length });
    } else {
        auto ops = read.cigar.retro.until!q{ !a.is_clipping };
        auto clipped = reduce!q{ a + b }(0, ops.map!q{ a.length });
        return read.position + read.basesCovered() + clipped;
    }
}

ushort computeScore(R)(auto ref R read) {
    return reduce!"a + b"(0, read.base_qualities.filter!"a >= 15").to!ushort;
}

auto collectSingleEndInfo(IndexedBamRead read, ReadGroupIndex read_group_index) {
    assert(read.ref_id != -1);
        
    SingleEndInfo result = void;
    result.coord = computeFivePrimeCoord(read);
    result.idx = read.index;
    result.score = computeScore(read);
    result.ref_id = cast(ushort)read.ref_id;
    result.reversed = read.is_reverse_strand ? 1 : 0;
    result.paired = (read.is_paired && !read.mate_is_unmapped) ? 1 : 0;

    auto rg = read_group_index.getId(getRG(read));
    result.library_id = cast(short)read_group_index.getLibraryId(rg);
    return result;
}

// may swap the two arguments
PairedEndsInfo combine(ref SingleEndInfo s1, ref SingleEndInfo s2) {
    assert(s1.library_id == s2.library_id);
    assert(s1.paired && s2.paired);
        
    if ((s2.ref_id < s1.ref_id) ||
        ((s2.ref_id == s1.ref_id) &&
         ((s2.coord < s1.coord) ||
          (s2.coord == s1.coord && s2.reversed < s1.reversed))))
        swap(s1, s2);

    PairedEndsInfo result = void;
    result.library_id = s1.library_id;
    result.ref_id1 = s1.ref_id;
    result.ref_id2 = s2.ref_id;
    result.reversed1 = s1.reversed;
    result.reversed2 = s2.reversed;
    result.coord1 = s1.coord;
    result.coord2 = s2.coord;
    result.score = cast(ushort)(s1.score + s2.score);
    result.idx1 = s1.idx;
    result.idx2 = s2.idx;
    return result;
}

// helpers for $(D samePosition)
SingleEndBasicInfo basicInfo(E)(auto ref E e) {
    static if (is(E == SingleEndBasicInfo))
        return e;
    else static if (is(E == SingleEndInfo))
             return e.basic_info;
    else static if (is(E == PairedEndsInfo))
             return e.read1_basic_info;
}
    
bool samePosition(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
    static if (is(E1 == PairedEndsInfo) && is(E2 == PairedEndsInfo)) {
        return *cast(ulong*)(&e1) == *cast(ulong*)(&e2) &&
            e1.ref_id2 == e2.ref_id2 && e1.coord2 == e2.coord2;
    } else {
        return basicInfo(e1).samePosition(basicInfo(e2));
    }
}

bool positionLessOrEq(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
    return !singleEndInfoComparator(basicInfo(e2), basicInfo(e1));
}
    
bool positionLess(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
    return !positionLessOrEq(e2, e1);
}

class SortedStorage(T, alias Cmp="a<b") {
    private {
        string _tmp_dir, _salt;
        string[] _filenames;
        MallocArray!T _buffer;
        TaskPool _pool;
        size_t _length;
    }

    this(string tmp_dir, TaskPool pool, size_t max_memory_usage) {
        _buffer = new MallocArray!T(max_memory_usage / T.sizeof);
        _tmp_dir = tmp_dir;
        _pool = pool;

        import std.file : dirEntries, SpanMode;
        _filenames = dirEntries(tmp_dir, T.stringof ~ "*",
                                SpanMode.shallow).map!(a => a.name).array;
        if (!_filenames.empty) { // this if-branch is for debug purposes only
            _salt = _filenames[0][T.stringof.length .. $][0 .. 4];
        } else {
            import std.random;
            char[4] tmp;
            auto gen = Random(unpredictableSeed);
            foreach (ref c; tmp[]) c = uniform!"[]"('a', 'z', gen);
            _salt = tmp[].idup;
        }
    }

    void put(T item) {
        if (_buffer.capacity == _buffer.data.length)
            flush();
        _buffer.put(item);
        ++_length;
    }

    size_t length() @property const {
        return _length;
    }

    void flush() {
        if (_buffer.data.length == 0)
            return;

        unstableSort!Cmp(_buffer.data, _pool);

        auto fn = buildPath(_tmp_dir, T.stringof ~ _salt ~
                            _filenames.length.to!string);

        auto compressor = new LZ4Compressor();
        auto lz4_output = std.stdio.File(fn, "w+");
        compressor.compress(cast(ubyte[])_buffer.data, lz4_output, 1);
        lz4_output.close();

        _filenames ~= fn;
        _buffer.reset();
    }

    void close(bool remove_files=false) {
        flush();
        _buffer.free();

        if (remove_files)
            removeTemporaryFiles();
    }

    void removeTemporaryFiles() {
        foreach (fn; _filenames)
            std.file.remove(fn);
    }

    auto reader() {
        SimpleReader!(T, LZ4File)[] readers;
        foreach (fn; _filenames)
            readers ~= new SimpleReader!(T, LZ4File)(fn);
        return nWayUnion!Cmp(readers);
    }
}

class SimpleReader(T, U) {
    this(string filename) {
        _file = U(filename);
        _buf = new ubyte[T.sizeof * 1024];
        popFront();
    }

    private {
        U _file;
        T _front;
        bool _empty;
        ubyte[] _buf;
        ubyte[] _raw_data;
    }

    bool empty() @property const { return _empty; }
    T front() @property const { return _front; }

    void popFront() {
        import std.algorithm : copy;

        while (_raw_data.length < T.sizeof) {
            auto remaining = _raw_data.length;
            copy(_raw_data[], _buf[0 .. remaining]);
            _raw_data = _file.rawRead(_buf[remaining .. $]);
            _raw_data = _buf[0 .. remaining + _raw_data.length];

            if (_raw_data.length == 0) {
                _empty = true;
                _file.close();
                return;
            }
        }
        _front = *(cast(T*)_raw_data.ptr);
        _raw_data = _raw_data[T.sizeof .. $];
    }
}

alias PEStorage = SortedStorage!(PairedEndsInfo, pairedEndsInfoComparator);
alias SEStorage = SortedStorage!(SingleEndInfo, singleEndInfoComparator);
alias SEBStorage = SortedStorage!(SingleEndBasicInfo, singleEndInfoComparator);

auto collectDuplicates(PEStorage pe_storage,
                       SEStorage se_storage,
                       SEBStorage pos_storage,
                       MarkDuplicatesConfig cfg,
                       TaskPool pool)
{
    auto pe_tmp = new MallocArray!PairedEndsInfo();
    auto se_tmp = new MallocArray!SingleEndInfo();

    auto pe = pe_storage.reader;
    auto se = se_storage.reader;
    auto pos = pos_storage.reader;

    auto duplicate_indices = new SortedStorage!ulong(cfg.tmpdir, pool, cfg.bufsize);

    while (true) {
        // sort of three-way merge

        // number of elements processed at the current position
        size_t pe_proc;
        size_t se_proc;

        pe_tmp.reset(); se_tmp.reset();

        // process paired ends
        if (!pe.empty && (se.empty || positionLessOrEq(pe.front, se.front))) {
            // this is essentially group-by operation,
            // but maximum score is computed on the fly for efficiency
            pe_tmp.put(pe.front);
            pe.popFront();
            size_t k = 1;
            size_t best_k = 0;
            while (!pe.empty) {
                if (!samePosition(pe.front, pe_tmp.data[0]))
                    break;
                if (pe.front.score > pe_tmp.data[best_k].score)
                    best_k = k;
                pe_tmp.put(pe.front);
                pe.popFront();
                ++k;
            }

            if (cfg.pe_callback !is null)
                cfg.pe_callback(pe_tmp.data);
            
            for (size_t i = 0; i < k; ++i) {
                if (i != best_k) {
                    auto paired_ends = pe_tmp.data[i];
                    duplicate_indices.put(paired_ends.idx1);
                    duplicate_indices.put(paired_ends.idx2);
                }
            }

            pe_proc = k;
        }

        // process single ends
        if (!se.empty && (pe_tmp.data.empty || positionLessOrEq(se.front, pe_tmp.data[0]))) {
            while (!pos.empty && positionLess(pos.front, se.front))
                pos.popFront();

            bool seen_paired_read = false;
            bool seen_fragment = false;

            auto cur = se.front;
            size_t k = 0;
            while (!se.empty) {
                if (!samePosition(se.front, cur))
                    break;
                if (se.front.paired)
                    seen_paired_read = true;
                else
                    seen_fragment = true;
                se_tmp.put(se.front); ++k;
                se.popFront();
            }

            size_t k_pe, k_pos;
            while (k_pe < pe_tmp.data.length && samePosition(pe_tmp.data[k_pe], cur)) {
                ++k_pe;
            }
            while (!pos.empty && samePosition(pos.front, cur)) {
                ++k_pos;
                pos.popFront();
            }

            if (k_pe + k_pos > 0)
                seen_paired_read = true;

            size_t total = k + k_pe + k_pos;

            if (cfg.se_callback !is null)
                cfg.se_callback(se_tmp.data, seen_paired_read);

            if (total < 2 || !seen_fragment) { /* do nothing */ }
            else if (seen_paired_read) {
                // there was a paired read at this position
                // => mark all single reads as duplicates
                for (size_t i = 0; i < k; ++i)
                    if (!se_tmp.data[i].paired)
                        duplicate_indices.put(se_tmp.data[i].idx);
            } else {
                size_t best_i = 0;
                for (size_t i = 0; i < k; ++i) {
                    if (se_tmp.data[i].score > se_tmp.data[best_i].score)
                        best_i = i;
                }
                for (size_t i = 0; i < k; ++i) {
                    if (i != best_i)
                        duplicate_indices.put(se_tmp.data[i].idx);
                }
            }
            se_proc = k;
        }

        if (se.empty && pe.empty) {
            break;
        }
            
        assert(pe_proc + se_proc > 0);
    }

    pe_tmp.free();
    se_tmp.free();

    duplicate_indices.close();
    return duplicate_indices;
}

auto getDuplicateOffsets(R)(R reads, ReadGroupIndex rg_index,
                            TaskPool pool, MarkDuplicatesConfig cfg) {

    auto single_ends = new SEStorage(cfg.tmpdir, pool, cfg.bufsize);
    auto paired_ends = new PEStorage(cfg.tmpdir, pool, cfg.bufsize);
    auto second_ends = new SEBStorage(cfg.tmpdir, pool, cfg.bufsize);

    scope(failure) {
        paired_ends.close(true);
        single_ends.close(true);
        second_ends.close(true);
    }

    size_t unmatched_pairs;

    foreach (pf; readPairsAndFragments(reads,
                                       cfg.hash_table_size_log2,
                                       cfg.overflow_list_size,
                                       cfg.tmpdir, pool)) {
        auto end1 = collectSingleEndInfo(pf.read1, rg_index);
        if (!pf.read2.isNull) {
            auto end2 = collectSingleEndInfo(pf.read2, rg_index);
            auto pair = combine(end1, end2);
            paired_ends.put(pair);
            second_ends.put(end2.basic_info);
        } else {
            if (end1.paired)
                ++unmatched_pairs;
            single_ends.put(end1);
        }
    }

    paired_ends.flush();
    single_ends.flush();
    second_ends.flush();

    stderr.writeln("  sorted ", paired_ends.length, " end pairs");

    stderr.writeln("     and ", single_ends.length, " single ends",
                   " (among them ", unmatched_pairs, " unmatched pairs)");

    StopWatch sw;
    stderr.write("  collecting indices of duplicate reads... ");
    sw.start();
    auto duplicates = collectDuplicates(paired_ends, single_ends, second_ends, cfg, pool);
    sw.stop(); stderr.writeln("  done in ", sw.peek().msecs, " ms"); sw.reset();
    stderr.writeln("  found ", duplicates.length, " duplicates");

    paired_ends.removeTemporaryFiles();
    single_ends.removeTemporaryFiles();
    second_ends.removeTemporaryFiles();

    return duplicates;
}

} /* end of algorithm */

void printUsage() {
    stderr.writeln("Usage: sambamba-markdup [options] <input.bam> [<input2.bam> [...]] <output.bam>");
    stderr.writeln("       By default, marks the duplicates without removing them");
    stderr.writeln();
    stderr.writeln("Options: -r, --remove-duplicates");
    stderr.writeln("                    remove duplicates instead of just marking them");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    number of threads to use");
    stderr.writeln("         -l, --compression-level=N");
    stderr.writeln("                    specify compression level of the resulting file (from 0 to 9)");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("                    show progressbar in STDERR");
    stderr.writeln("         --tmpdir=TMPDIR");
    stderr.writeln("                    specify directory for temporary files");
    stderr.writeln();
    stderr.writeln("Performance tweaking parameters");
    stderr.writeln("         --hash-table-size=HASH_TABLE_SIZE");
    stderr.writeln("                    size of hash table for finding read pairs (default is 262144 reads);");
    stderr.writeln("                    will be rounded down to the nearest power of two;");
    stderr.writeln("                    should be > (average coverage) * (insert size) for good performance");
    stderr.writeln("         --overflow-list-size=OVERFLOW_LIST_SIZE");
    stderr.writeln("                    size of the overflow list where reads, thrown from the hash table,");
    stderr.writeln("                    get a second chance to meet their pairs (default is 200000 reads);");
    stderr.writeln("                    increasing the size reduces the number of temporary files created");
    stderr.writeln("         --io-buffer-size=BUFFER_SIZE");
    stderr.writeln("                    two buffers of BUFFER_SIZE *megabytes* each are used");
    stderr.writeln("                    for reading and writing BAM during the second pass (default is 128)");
}

version(standalone) {
    int main(string[] args) {
        return markdup_main(args);
    }
}

int markdup_main(string[] args) {

    MarkDuplicatesConfig cfg;
    cfg.tmpdir = defaultTmpDir();

    bool remove_duplicates;
    size_t n_threads = totalCPUs;
    bool show_progress;
    size_t io_buffer_size = 128;
    size_t hash_table_size;
    int compression_level = -1;

    bool cmp_with_picard_mode; // for development purposes!
    
    StopWatch sw;
    sw.start();  
    
    try {
        getopt(args,
           std.getopt.config.caseSensitive,
           "remove-duplicates|r", &remove_duplicates,
           "nthreads|t", &n_threads,
           "compression-level|l", &compression_level,
           "show-progress|p", &show_progress,
           "tmpdir", &cfg.tmpdir,
           "hash-table-size", &hash_table_size,
           "overflow-list-size", &cfg.overflow_list_size,
           "io-buffer-size", &io_buffer_size,
           "compare-with-picard-mode", &cmp_with_picard_mode);

        if (args.length < 3) {
            printUsage();
            return 0;
        }

        foreach (arg; args[1 .. $-1])
            protectFromOverwrite(arg, args[$-1]);
        cfg.tmpdir = randomSubdir(cfg.tmpdir, "markdup-");

        auto pool = new TaskPool(n_threads);
        scope(exit) pool.finish();

        if (cmp_with_picard_mode) {
            static class PicardChecker {
                static string output_filename;
                this(MarkDuplicatesConfig cfg, string fn1, string fn2, TaskPool pool) {
                    output_filename = buildPath(cfg.tmpdir, "_diff1.bam");
                    auto output_filename2 = buildPath(cfg.tmpdir, "_diff2.bam");
                    auto b1 = new BamReader(fn1, pool);
                    auto b2 = new BamReader(fn2, pool);
                    b1.setBufferSize(64_000_000);
                    b2.setBufferSize(64_000_000);
                    auto w = new BamWriter(output_filename, 1, pool);
                    auto w2 = new BamWriter(output_filename2, 1, pool);
                    scope(exit) {
                        w.finish();
                        w2.finish();
                        stderr.writeln("Saved differing reads to ", output_filename,
                                       " and ", output_filename2);
                    }
                    w.writeSamHeader(b1.header);
                    w.writeReferenceSequenceInfo(b1.reference_sequences);
                    w2.writeSamHeader(b2.header);
                    w2.writeReferenceSequenceInfo(b2.reference_sequences);
                    foreach (pair; zip(b1.reads, b2.reads)) {
                        if (pair[0].is_duplicate != pair[1].is_duplicate) {
                            w.writeRecord(pair[0]);
                            w2.writeRecord(pair[1]);
                        }
                    }
                }

                void check(PairedEndsInfo[] dups) {
                    if (dups.length != 2 || dups[0].score != dups[1].score) {
                        writefln("weird group of PE duplicates found, their indices: %(%s, %); scores: %(%s, %)",
                                 roundRobin(dups.map!"a.idx1", dups.map!"a.idx2"),
                                 dups.map!"a.score");
                    }
                }

                void check(SingleEndInfo[] dups, bool has_paired) {
                    if (has_paired || dups.length != 2 || dups[0].score != dups[1].score) {
                        writefln("weird group of SE duplicates found, their indices: %(%s, %); scores: %(%s, %)",
                                 dups.map!"a.idx", dups.map!"a.score");
                    }
                }
            }
            auto checker = new PicardChecker(cfg, args[1], args[2], pool);
            args[1] = checker.output_filename;
            args[2] = "/dev/null";
            cfg.pe_callback = (pe_dups) => checker.check(pe_dups);
            cfg.se_callback = (se_dups, has_paired) => checker.check(se_dups, has_paired);
        }

        io_buffer_size <<= 20; // -> convert to megabytes

        cfg.hash_table_size_log2 = 10; // FIXME: overrides default value of 18
        while ((2UL << cfg.hash_table_size_log2) <= hash_table_size)
            cfg.hash_table_size_log2 += 1;
        // 2^^(cfg.hash_table_size_log2 + 1) > hash_table_size

        // Set up the BAM reader and pass in the thread pool
        auto bam = new MultiBamReader(args[1 .. $-1], pool);
        auto n_refs = bam.reference_sequences.length;
        enforce(n_refs < 16384, "More than 16383 reference sequences are unsupported");

        auto rg_index = new ReadGroupIndex(bam.header);

        stderr.writeln("finding positions of the duplicate reads in the file...");

        InputRange!IndexedBamRead reads;
        shared(ProgressBar) bar;

        void initInputs() {
            if (!show_progress)
                reads = bam.reads.withIndices.inputRangeObject;
            else {
                bar = new shared(ProgressBar)();
                reads = bam.readsWithProgress((lazy float p) { bar.update(p); },
                                              () { bar.finish(); }).withIndices
                           .inputRangeObject;
            }
        }

        initInputs();
        auto dup_idx_storage = getDuplicateOffsets(reads, rg_index, pool, cfg);

        auto elapsed = sw.peek();
        stderr.writeln("collected list of positions in ",
                       elapsed.seconds / 60, " min ",
                       elapsed.seconds % 60, " sec");

        // marking or removing duplicates
        bam = new MultiBamReader(args[1 .. $-1], pool);
        bam.setBufferSize(io_buffer_size);
        auto out_stream = new BufferedFile(args[$-1], FileMode.OutNew, io_buffer_size);
        auto writer = new BamWriter(out_stream, compression_level, pool);
        writer.setFilename(args[$-1]);
        scope(exit) writer.finish();
        writer.writeSamHeader(bam.header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);

        stderr.writeln(remove_duplicates ? "removing" : "marking", " duplicates...");
        
        initInputs();

        auto indices = dup_idx_storage.reader;

        foreach (read; reads) {
            if (!indices.empty && read.index == indices.front) {
                assumeUnique(read).is_duplicate = true;
                indices.popFront();
            } else if (!read.is_secondary_alignment && !read.is_supplementary) {
                assumeUnique(read).is_duplicate = false;
            }

            if (read.is_duplicate && remove_duplicates)
                continue;

            writer.writeRecord(read);
        }

        dup_idx_storage.removeTemporaryFiles();
        try {
          std.file.rmdirRecurse(cfg.tmpdir);
        } catch (FileException e) {
        }

        sw.stop();
        stderr.writeln("total time elapsed: ",
                       sw.peek().seconds / 60, " min ",
                       sw.peek().seconds % 60, " sec");

    } catch (Throwable e) {
        stderr.writeln("sambamba-markdup: ", e.msg);
        version(development) { throw e; }
        return 1;
    }

    return 0;
}
