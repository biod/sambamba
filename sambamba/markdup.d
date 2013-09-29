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
module sambamba.markdup;

import std.stdio;
import std.getopt;
import sambamba.utils.common.progressbar;
import thirdparty.unstablesort;

import bio.bam.reader, bio.bam.readrange, bio.bam.writer, bio.bam.referenceinfo,
       bio.bam.read, bio.sam.header, bio.bam.abstractreader,
       bio.core.bgzf.virtualoffset;
import std.traits, std.typecons, std.range, std.algorithm, std.parallelism, 
       std.exception, std.file, std.typetuple, std.conv, std.array, std.bitmanip,
       std.c.stdlib, std.datetime, std.stream : BufferedFile, FileMode;

struct ReadPair {
    BamReadBlock read1;
    BamReadBlock read2;
    this(BamReadBlock r1, BamReadBlock r2) {
        read1 = r1;
        read2 = r2;
    }
}

struct ReadPairOrFragment {
    BamReadBlock read1;
    Nullable!BamReadBlock read2;

    this(BamReadBlock r1) {
        read1 = r1;
    }

    this(BamReadBlock r1, BamReadBlock r2) {
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
}

/// 48 bytes; fast access to read group;
/// also stores precomputed hash.
struct HReadBlock {
    BamReadBlock read;
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
static assert(HReadBlock.sizeof == 48);

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
auto _rse2brb(T)(T t) { return BamReadBlock(t[1], t[2], t[0]); }
auto _rse2rs(T)(T t) { return map!_rse2brb(zip(t[0].reads, t[1], t[2])); }
auto _rse2r(BamReader[] r, VirtualOffset[][] s, VirtualOffset[][] e) {
    return map!_rse2rs(zip(r, s, e)).array()
              .nWayUnion!compareReadNamesAndReadGroups();
}

private auto readsFromTempFiles(size_t buf_size, string[] tmp_filenames,
                                TaskPool pool) {
    BamReader[] readers;
    VirtualOffset[][] start_virtual_offsets;
    VirtualOffset[][] end_virtual_offsets;
    foreach (fn; tmp_filenames) {
        readers ~= new BamReader(fn, pool);
        auto offsets = cast(VirtualOffset[])std.file.read(fn ~ ".vo");
        start_virtual_offsets ~= offsets[0 .. $ / 2];
        end_virtual_offsets ~= offsets[$ / 2 .. $];
        readers[$ - 1].setBufferSize(buf_size / tmp_filenames.length);
    }

    return _rse2r(readers, start_virtual_offsets, end_virtual_offsets);
}

struct CollateReadPairRange(R, bool keepFragments, alias charsHashFunc)
    if (isInputRange!R && is(Unqual!(ElementType!R) == BamReadBlock))
{
    private {
        static auto wrapper(R reads, TaskPool pool) {
            auto r1 = reads.until!q{ a.ref_id == -1 }
                           .filter!q{ !a.is_unmapped }
                           .filter!q{ !a.is_secondary_alignment };

            static if (keepFragments) {
                auto r2 = r1;
            } else {
                auto r2 = r1.filter!q{a.is_paired && !a.mate_is_unmapped}();
            }

            return pool.map!(makeHReadBlock!charsHashFunc)(r2, 1024);
        }

        ReturnType!wrapper _reads;
        TaskPool _task_pool;

        alias BamReadBlock Read;

        HReadBlock[] _table; // size is always power of two
        size_t _table_mask;

        ulong _min_vo, _max_vo;

        version(profile) {
            size_t _duped_during_compaction;
            StopWatch _compact_sw;
        }

        void _compact() {
            enum max_diff = 10_000_000UL << 16;
            if (_max_vo - _min_vo < max_diff * 6)
                return;

            version(profile) { _compact_sw.start(); scope(exit) _compact_sw.stop(); }
            
            _min_vo = _max_vo;
            foreach (ref r; _table)
                if (!r.isNull && r.is_slice_backed
                    && _max_vo - cast(ulong)r.start_virtual_offset > max_diff)
                {
                    version(profile) ++_duped_during_compaction;
                    r = r.dup;
                } else if (!r.isNull && r.is_slice_backed) {
                    _min_vo = min(_min_vo, cast(ulong)r.start_virtual_offset);
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
        VirtualOffset[] _tmp_s_vo;
        VirtualOffset[] _tmp_e_vo;
        Nullable!Read _tmp_r1;
        ReturnType!readsFromTempFiles _tmp_reads;
    }

    this(R reads, ubyte table_size_log2, size_t overflow_list_size,
         string tmp_dir = "/tmp", TaskPool task_pool = taskPool)
    {
        enforce(overflow_list_size > 0);
        _tmp_dir = tmp_dir;

        _reads = wrapper(reads, task_pool);
        _task_pool = task_pool;
        setSource(Source.hashTable);

        _table = new HReadBlock[1 << table_size_log2];
        _table_mask = (1 << table_size_log2) - 1;
        _overflow_list = new HReadBlock[overflow_list_size];

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
        _max_vo = max(_max_vo, cast(ulong)read.start_virtual_offset);
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
        _tmp_w.writeSamHeader(_reader.header);
        _tmp_w.writeReferenceSequenceInfo(_reader.reference_sequences);
        _tmp_written = 0;
    }

    void closeTmpWriter() {
        if (_tmp_w !is null)
            _tmp_w.finish();
        auto vo_fn = _tmp_filenames[$ - 1] ~ ".vo";
        std.file.write(vo_fn, cast(ubyte[])_tmp_s_vo);
        std.file.append(vo_fn, cast(ubyte[])_tmp_e_vo);
        _tmp_s_vo.length = 0;
        _tmp_e_vo.length = 0;
        assumeSafeAppend(_tmp_s_vo);
        assumeSafeAppend(_tmp_e_vo);
    }

    void dumpTmpRecord(R)(auto ref R read) {
        assert(_tmp_w !is null);
        _tmp_w.writeRecord(read);
        _tmp_written++;
        _tmp_s_vo ~= read.start_virtual_offset;
        _tmp_e_vo ~= read.end_virtual_offset;
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
            std.file.remove(fn ~ ".vo");
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
(R reads, ubyte table_size_log2=18, size_t overflow_list_size=200000,
 string tmp_dir="/tmp", TaskPool task_pool=taskPool) {
    return CollateReadPairRange!(R, false, hashFunc)
            (reads, table_size_log2, overflow_list_size, tmp_dir, task_pool);
}

auto readPairsAndFragments(alias hashFunc=simpleHash, R)
(R reads, ubyte table_size_log2=18, size_t overflow_list_size=200000,
 string tmp_dir="/tmp", TaskPool task_pool=taskPool) {
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
    VirtualOffset vo;
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
    VirtualOffset vo1, vo2;

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

// Assumptions needed to save memory.
// Trivial ones:
static assert(2 * VirtualOffset.sizeof <= PairedEndsInfo.sizeof);
static assert(VirtualOffset.sizeof <= SingleEndInfo.sizeof);
// Not so trivial:
//
// Suppose we have A pairs and B single reads.
// Our goal is to place virtual offsets of duplicates into one of
// already existing chunks of memory. The two already existing arrays
// consume A * P and B * S bytes where P and S are the sizeofs.
// In order to reuse the memory, we need
// (2A + B) * V <= max(A * P, B * S)  (where V = VirtualOffset.sizeof)
// Suppose the contrary: (2A + B) * V > A * P and (2A + B) * V > B * S.
// This implies B * V > A * (P - 2V), A * 2V > B * (S - V), or
// V / (P - 2V) > A/B > (S - V) / 2V (in real arithmetic).
// Therefore, if V / (P - 2V) <= (S - V) / 2V, everything is fine.
// This is equivalent to 2V^2 <= (S - V) * (P - 2V)
//                       (S - V) * P >= 2 * S * V
static assert(SingleEndInfo.sizeof > VirtualOffset.sizeof);
static assert((SingleEndInfo.sizeof - VirtualOffset.sizeof)
                                    * PairedEndsInfo.sizeof
              >= 2 * SingleEndInfo.sizeof * VirtualOffset.sizeof);

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
    string tmpdir = "/tmp";
}

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

/// The returned result must be freed using std.c.malloc.free
VirtualOffset[] getDuplicateOffsets(R)(R reads, ReadGroupIndex rg_index,
                                       TaskPool pool, MarkDuplicatesConfig cfg) {

    static int computeFivePrimeCoord(R)(auto ref R read) {
        if (!read.is_reverse_strand) {
            auto ops = read.cigar.until!q{ !a.is_clipping };
            return read.position - reduce!q{ a + b }(0, ops.map!q{ a.length });
        } else {
            auto ops = read.cigar.retro.until!q{ !a.is_clipping };
            auto clipped = reduce!q{ a + b }(0, ops.map!q{ a.length });
            return read.position + read.basesCovered() + clipped;
        }
    }

    static ushort computeScore(R)(auto ref R read) {
        return reduce!"a + b"(0, read.base_qualities.filter!"a >= 15").to!ushort;
    }

    static auto collectSingleEndInfo(BamReadBlock read,
                                     ReadGroupIndex read_group_index) {
        assert(read.ref_id != -1);
        
        SingleEndInfo result = void;
        result.coord = computeFivePrimeCoord(read);
        result.vo = read.start_virtual_offset;
        result.score = computeScore(read);
        result.ref_id = cast(ushort)read.ref_id;
        result.reversed = read.is_reverse_strand ? 1 : 0;
        result.paired = (read.is_paired && !read.mate_is_unmapped) ? 1 : 0;

        auto rg = read_group_index.getId(getRG(read));
        result.library_id = cast(short)read_group_index.getLibraryId(rg);
        return result;
    }

    // may swap the two arguments
    static PairedEndsInfo combine(ref SingleEndInfo s1,
                                  ref SingleEndInfo s2)
    {
        assert(s1.library_id == s2.library_id);
        assert(s1.paired && s2.paired);
        
        if (s2.ref_id <= s1.ref_id)
            if (s2.ref_id < s1.ref_id || s2.coord < s1.coord)
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
        result.vo1 = s1.vo;
        result.vo2 = s2.vo;
        return result;
    }

    // helpers for $(D samePosition)
    static SingleEndBasicInfo basicInfo(E)(auto ref E e) {
        static if (is(E == SingleEndBasicInfo))
            return e;
        else static if (is(E == SingleEndInfo))
            return e.basic_info;
        else static if (is(E == PairedEndsInfo))
            return e.read1_basic_info;
    }
    
    static bool samePosition(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
        static if (is(E1 == PairedEndsInfo) && is(E2 == PairedEndsInfo)) {
            return *cast(ulong*)(&e1) == *cast(ulong*)(&e2) &&
                e1.ref_id2 == e2.ref_id2 && e1.coord2 == e2.coord2;
        } else {
            return basicInfo(e1).samePosition(basicInfo(e2));
        }
    }

    static bool positionLessOrEq(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
        return !singleEndInfoComparator(basicInfo(e2), basicInfo(e1));
    }
    
    static bool positionLess(E1, E2)(auto ref E1 e1, auto ref E2 e2) {
        return !positionLessOrEq(e2, e1);
    }

    /// Frees one of $(D _pe) or $(D _se). Reuses memory of the other.
    /// Returns an array pointing to data which is to be freed after use.
    static VirtualOffset[] collectDuplicates(MallocArray!PairedEndsInfo _pe,
                                             MallocArray!SingleEndInfo _se,
                                             MallocArray!SingleEndBasicInfo _pos)
    {
        // yeah, that's dirty, but allows to minimize memory usage
        auto pvo = cast(VirtualOffset*)cast(void*)_pe.data.ptr;
        auto svo = cast(VirtualOffset*)cast(void*)_se.data.ptr;
        size_t n_vo_pe, n_vo_se; // how much VOs we store in each array

        auto pe = _pe.data;
        auto se = _se.data;
        auto pos = _pos.data;

        size_t pe_total_mem = pe.length * PairedEndsInfo.sizeof;
        size_t se_total_mem = se.length * SingleEndInfo.sizeof;

        // sort of three-way merge
        while (true) {
            // number of elements processed at the current position
            size_t pe_proc;
            size_t se_proc;

            // process paired ends
            if (!pe.empty && (se.empty || positionLessOrEq(pe.front, se.front))) {
                // this is essentially group-by operation,
                // but maximum score is computed on the fly for efficiency
                auto cur = pe.front;
                size_t k = 0;
                size_t best_k = 0;
                ++k;
                while (k < pe.length && samePosition(pe[k], cur)) {
                    if (pe[k].score > pe[best_k].score)
                        best_k = k;
                    ++k;
                }
                // virtual offsets are places to the chunk of memory
                // that we already processed => it's safe because
                // 2 * VirtualOffset.sizeof <= PairedEndsInfo.sizeof
                for (size_t i = 0; i < k; ++i) {
                    if (i != best_k) {
                        pvo[n_vo_pe++] = pe[i].vo1;
                        pvo[n_vo_pe++] = pe[i].vo2;
                    }
                }

                pe_proc = k;
            }

            // process single ends
            if (!se.empty && (pe.empty || positionLessOrEq(se.front, pe.front))) {
                while (!pos.empty && positionLess(pos.front, se.front))
                    pos.popFront();

                bool seen_paired_read = false;
                bool seen_fragment = false;

                auto cur = se.front;
                size_t k = 0;
                while (k < se.length && samePosition(se[k], cur)) {
                    if (se[k].paired)
                        seen_paired_read = true;
                    else
                        seen_fragment = true;
                    ++k;
                }

                size_t k_pe, k_pos;
                while (k_pe < pe.length && samePosition(pe[k_pe], cur))
                    ++k_pe;
                while (k_pos < pos.length && samePosition(pos[k_pos], cur))
                    ++k_pos;

                if (k_pe + k_pos > 0)
                    seen_paired_read = true;

                size_t total = k + k_pe + k_pos;

                if (total < 2 || !seen_fragment) { /* do nothing */ }
                else if (seen_paired_read) {
                    // there was a paired read at this position
                    // => mark all single reads as duplicates
                    for (size_t i = 0; i < k; ++i)
                        if (!se[i].paired)
                            svo[n_vo_se++] = se[i].vo;
                } else {
                    size_t best_i = 0;
                    for (size_t i = 0; i < k; ++i) {
                        if (se[i].score > se[best_i].score)
                            best_i = i;
                        ++i;
                    }
                    for (size_t i = 0; i < k; ++i) {
                        if (i != best_i)
                            svo[n_vo_se++] = se[i].vo;
                    }
                }
                se_proc = k;
            }

            // move to the next position
            pe = pe[pe_proc .. $];
            se = se[se_proc .. $];

            if (se.empty && pe.empty)
                break;
            
            assert(pe_proc + se_proc > 0);
        }

        // finished collecting duplicates but they are in two different arrays;
        // now we can join them.
        size_t pe_used = n_vo_pe * VirtualOffset.sizeof;
        size_t se_used = n_vo_se * VirtualOffset.sizeof;

        VirtualOffset* p;
        size_t n = n_vo_pe + n_vo_se;
        if (pe_used + se_used <= pe_total_mem) {
            pvo[n_vo_pe .. n] = svo[0 .. n_vo_se];
            _se.free();
            p = pvo;
        } else {
            svo[n_vo_se .. n_vo_pe + n_vo_se] = pvo[0 .. n_vo_pe];
            _pe.free();
            p = svo;
        }
        std.c.stdlib.realloc(p, n * VirtualOffset.sizeof);
        return p[0 .. n];
    }

    auto single_ends = new MallocArray!SingleEndInfo();
    auto paired_ends = new MallocArray!PairedEndsInfo();
    auto second_ends = new MallocArray!SingleEndBasicInfo();

    scope(failure) {
        paired_ends.free();
        single_ends.free();
    }

    scope(exit) second_ends.free();

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
            single_ends.put(end1);
        }
    }
    
    StopWatch sw;

    stderr.write("  sorting ", paired_ends.data.length, " end pairs... ");
    sw.start();
    unstableSort!pairedEndsInfoComparator(paired_ends.data, true);
    unstableSort!singleEndInfoComparator(second_ends.data, true);
    sw.stop(); stderr.writeln("  done in ", sw.peek().msecs, " ms"); sw.reset();

    stderr.write("  sorting ", single_ends.data.length, " single ends",
                 " (among them ",
                 single_ends.data.filter!q{ a.paired }.walkLength(),
                 " unmatched pairs)... ");
    sw.start();
    unstableSort!singleEndInfoComparator(single_ends.data, true);
    sw.stop(); stderr.writeln("done in ", sw.peek().msecs, " ms"); sw.reset();

    stderr.write("  collecting virtual offsets of duplicate reads... ");
    sw.start();
    auto duplicates = collectDuplicates(paired_ends, single_ends, second_ends);
    sw.stop(); stderr.writeln("  done in ", sw.peek().msecs, " ms"); sw.reset();

    stderr.write("  found ", duplicates.length, " duplicates, sorting the list... ");
    sw.start();
    unstableSort(duplicates);
    sw.stop(); stderr.writeln("  done in ", sw.peek().msecs, " ms"); sw.reset();

    return duplicates;
}

void printUsage() {
    stderr.writeln("Usage: sambamba-markdup [options] <input.bam>");
    stderr.writeln("       By default, marks the duplicates without removing them");
    stderr.writeln();
    stderr.writeln("Options: -o, --output-filename=FILENAME");
    stderr.writeln("                    the name of the output file; required option");
    stderr.writeln("         -r, --remove-duplicates");
    stderr.writeln("                    remove duplicates instead of just marking them");
    stderr.writeln("         -t, --nthreads=NTHREADS");
    stderr.writeln("                    number of threads to use");
    stderr.writeln("         -l, --compression-level=N");
    stderr.writeln("                    specify compression level of the resulting file (from 0 to 9)");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("                    show progressbar in STDERR");
    stderr.writeln("         --tmpdir=TMPDIR");
    stderr.writeln("                    specify directory for temporary files; default is /tmp");
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
    string output_filename;
    bool remove_duplicates;
    size_t n_threads = totalCPUs;
    bool show_progress;
    size_t io_buffer_size = 128;
    size_t hash_table_size;
    int compression_level = -1;
    
    getopt(args,
           std.getopt.config.caseSensitive,
           "output-filename|o", &output_filename,
           "remove-duplicates|r", &remove_duplicates,
           "nthreads|t", &n_threads,
           "compression-level|l", &compression_level,
           "show-progress|p", &show_progress,
           "tmpdir", &cfg.tmpdir,
           "hash-table-size", &hash_table_size,
           "overflow-list-size", &cfg.overflow_list_size,
           "io-buffer-size", &io_buffer_size);

    if (args.length == 1) {
        printUsage();
        return 0;
    }

    if (output_filename is null) {
        stderr.writeln("Please specify output filename");
        return 1;
    }

    io_buffer_size <<= 20; // -> megabytes

    cfg.hash_table_size_log2 = 10;
    while ((2 << cfg.hash_table_size_log2) < hash_table_size)
        cfg.hash_table_size_log2 += 1;

    try {
        auto pool = new TaskPool(n_threads);
        scope(exit) pool.finish();

        auto bam = new BamReader(args[1], pool);
        auto n_refs = bam.reference_sequences.length;
        enforce(n_refs < 16384, "More than 16383 reference sequences are unsupported");

        auto rg_index = new ReadGroupIndex(bam.header);

        VirtualOffset[] offsets;
        scope(exit) std.c.stdlib.free(offsets.ptr);

        StopWatch sw;
        stderr.writeln("finding positions of the duplicate reads in the file...");
        sw.start();
        if (!show_progress)
            offsets = getDuplicateOffsets(bam.reads!withOffsets(), rg_index, pool, cfg);
        else {
            auto bar = new shared(ProgressBar)();
            auto reads = bam.readsWithProgress!withOffsets((lazy float p) { bar.update(p); },
                                                           () { bar.finish(); });
            offsets = getDuplicateOffsets(reads, rg_index, pool, cfg);
        }
        sw.stop();
        stderr.writeln("collected list of positions in ",
                       sw.peek().seconds / 60, " min ",
                       sw.peek().seconds % 60, " sec");

        // marking or removing duplicates
        bam = new BamReader(args[1], pool);
        bam.setBufferSize(io_buffer_size);
        auto out_stream = new BufferedFile(output_filename, FileMode.OutNew, io_buffer_size);
        auto writer = new BamWriter(out_stream, compression_level, pool);
        scope(exit) writer.finish();
        writer.writeSamHeader(bam.header);
        writer.writeReferenceSequenceInfo(bam.reference_sequences);

        stderr.writeln(remove_duplicates ? "removing" : "marking", " duplicates...");
        sw.start();
        
        InputRange!BamReadBlock reads;
        if (!show_progress) {
            reads = inputRangeObject(bam.reads!withOffsets());
        } else {
            auto bar = new shared(ProgressBar)();
            reads = inputRangeObject(bam.readsWithProgress!withOffsets((lazy float p) { bar.update(p); },
                                                                       () { bar.finish(); }));
        }

        size_t k = 0;
        foreach (read; reads) {
            if (k < offsets.length && read.start_virtual_offset == offsets[k])
                {
                    ++k;
                    assumeUnique(read).is_duplicate = true;
                    if (remove_duplicates)
                        continue;
                } else {
                assumeUnique(read).is_duplicate = false;
            }
            writer.writeRecord(read);
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
