import bamfile;
import samheader;
import reference;
import alignment;
import tagvalue;
import validation.alignment;

import core.memory : GC;
import core.runtime : Runtime;
import std.c.stdlib;
import std.string;
import std.range;
import std.algorithm;
import core.stdc.string : memcpy;
import std.conv;

import std.traits;

extern(C) void libbam_init() {
    Runtime.initialize();
    GC.disable();
}

void main() {}

static string last_error;

extern(C) immutable(char)* 
get_last_error() {
    return last_error.ptr;
}

extern(C) void* 
bamfile_new(const char* filename) {
    try {
        BamFile* f = cast(BamFile*)malloc(BamFile.sizeof);
        scope(failure) free(f);

        emplace(f, to!string(filename));
        GC.addRange(f, BamFile.sizeof);
        return f;
    } catch (Throwable e) {
        last_error = e.msg;
        return null;
    }
}

extern(C) bool
bamfile_has_index(BamFile* f) {
    return f.has_index;
}

extern(C) void 
bamfile_destroy(BamFile* f) {
    f.close();
    GC.removeRange(f);
    free(cast(void*)f);
}

extern(C) void
bamfile_rewind(BamFile* f) {
    f.rewind();
}

extern(C) SamHeader 
bamfile_get_header(BamFile* f) {
    return f.header;
}

struct Array(T) {
    size_t length;
    T* ptr;
}

extern(C) Array!ReferenceSequenceInfo 
bamfile_get_reference_sequences(BamFile* f) {
    typeof(return) arr;
    arr.length = f.reference_sequences.length;
    arr.ptr = f.reference_sequences.ptr;
    return arr;
}

alias InputRange!Alignment AlignmentRange;

extern(C) AlignmentRange* 
bamfile_get_alignments(BamFile* f) {
    AlignmentRange* ar = cast(AlignmentRange*)malloc(AlignmentRange.sizeof);
    AlignmentRange _ar = inputRangeObject(f.alignments);
    memcpy(ar, &_ar, _ar.sizeof);
    GC.addRange(ar, AlignmentRange.sizeof);
    return ar;
}

extern(C) AlignmentRange* 
bamfile_get_valid_alignments(BamFile* f) {
    AlignmentRange* ar = cast(AlignmentRange*)malloc(AlignmentRange.sizeof);
    auto range = filter!((Alignment a){return isValid(a);})(f.alignments);
    AlignmentRange _ar = inputRangeObject(range);
    memcpy(ar, &_ar, _ar.sizeof);
    GC.addRange(ar, AlignmentRange.sizeof);
    return ar;
}

extern(C) AlignmentRange*
bamfile_fetch_alignments(BamFile* f, const char* chr, int beg, int end) {
    AlignmentRange* ar = cast(AlignmentRange*)malloc(AlignmentRange.sizeof);
    auto range = (*f)[to!string(chr)][beg .. end];
    AlignmentRange _ar = inputRangeObject(range);
    memcpy(ar, &_ar, _ar.sizeof);
    GC.addRange(ar, AlignmentRange.sizeof);
    return ar;
}

extern(C) void 
alignment_range_destroy(AlignmentRange* ar) {
    GC.removeRange(ar);
    free(cast(void*)ar);
}

extern(C) Alignment* 
alignment_range_front(AlignmentRange* ar) {
    Alignment* a = cast(Alignment*)malloc(Alignment.sizeof);
    auto _a = ar.front;
    memcpy(a, &_a, _a.sizeof);
    GC.addRange(ar, Alignment.sizeof);
    return a;
}

/// returns 1 if range is empty after popfront,
///         0 if non-empty,
///        -1 if error occurred
extern(C) int
alignment_range_popfront(AlignmentRange* ar) {
    try {
        ar.popFront();
        return ar.empty() ? 1 : 0;
    } catch (Throwable e) {
        last_error = e.msg;
        return -1;
    }
}

extern(C) bool 
alignment_range_empty(AlignmentRange* ar) {
    return ar.empty();
}

extern(C) void 
alignment_destroy(Alignment* a) {
    GC.removeRange(a);
    free(cast(void*)a);
}

extern(C) bool
alignment_is_valid(Alignment* a) {
    return isValid(*a);
}

extern(C) int 
alignment_ref_id(Alignment* a) {
    return a.ref_id;
}

extern(C) int 
alignment_position(Alignment* a) {
    return a.position;
}

extern(C) ushort 
alignment_bin(Alignment* a) {
    return a.bin;
}

extern(C) ubyte 
alignment_mapping_quality(Alignment* a) {
    return a.mapping_quality;
}

extern(C) ushort 
alignment_flag(Alignment* a) {
    return a.flag;
}

extern(C) int 
alignment_sequence_length(Alignment* a) {
    return a.sequence_length;
}

extern(C) int 
alignment_next_ref_id(Alignment* a) {
    return a.next_ref_id;
}

extern(C) int 
alignment_next_pos(Alignment* a) {
    return a.next_pos;
}

extern(C) int 
alignment_template_length(Alignment* a) {
    return a.template_length;
}

extern(C) Array!(immutable(char))
alignment_read_name(Alignment* a) {
    typeof(return) arr;
    arr.length = a.read_name.length;
    arr.ptr = a.read_name.ptr;
    return arr;
}

extern(C) immutable(char)*
alignment_cigar_string(Alignment* a) {
    return toStringz(a.cigar_string());
}

extern(C) immutable(char)*
alignment_sequence(Alignment* a) {
    return toStringz(to!string(a.sequence()));
}

extern(C) Array!ubyte
alignment_phred_base_quality(Alignment* a) {
    Array!ubyte arr;
    arr.length = a.phred_base_quality.length;
    arr.ptr = cast(ubyte*)(a.phred_base_quality.ptr);
    return arr;
}

extern(C) Value*
alignment_get_tag_value(Alignment* a, char* str) {
    Value* v = cast(Value*)malloc(Value.sizeof);
    auto _v = a.tags[to!string(str)];
    memcpy(v, &_v, _v.sizeof);
    GC.addRange(v, Value.sizeof);
    return v;
}

extern(C) void
tag_value_destroy(Value* v) {
    GC.removeRange(v);
    free(v);
}

struct DHash {
    string[] keys = void;
    Value[] values = void;
}

extern(C) void
dhash_destroy(DHash* hash) {
    GC.removeRange(hash);
    free(hash);
}

extern(C) DHash*
alignment_get_all_tags(Alignment* a) {
    DHash* hash = cast(DHash*)malloc(DHash.sizeof);
    GC.addRange(hash, DHash.sizeof);
    hash.keys = new string[0];
    hash.values = new Value[0];
    foreach (k, v; a.tags) {
        hash.keys ~= k;
        hash.values ~= v;
    }
    return hash;
}
