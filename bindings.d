import bamfile;
import samheader;

import core.memory : GC;
import core.runtime : Runtime;
import std.c.stdlib;
import std.string;
import core.stdc.string : memcpy;
import std.conv;

extern(C) void libbam_init() {
	Runtime.initialize();
}

void main() {}

static string last_error;

extern(C) immutable(char)* get_last_error() {
	return last_error.ptr;
}

extern(C) void* bamfile_new(const char* filename) {
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

extern(C) void bamfile_destroy(BamFile* f) {
    f.close();
    GC.removeRange(f);
    free(cast(void*)f);
}

extern(C) SamHeader bamfile_get_header(BamFile* f) {
	return f.header;
}
