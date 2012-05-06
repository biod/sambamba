import bamfile;
import samheader;

import core.memory : GC;
import core.runtime : Runtime;
import std.c.stdlib;
import core.stdc.string : memcpy;
import std.conv;

extern(C) void libbam_init() {
	Runtime.initialize();
}

void main() {}

extern(C) void* bamfile_new(const char* filename) {
    BamFile* f = cast(BamFile*)malloc(BamFile.sizeof);
    emplace(f, to!string(filename));
    GC.addRange(f, BamFile.sizeof);
    return f;
}

extern(C) void bamfile_destroy(BamFile* f) {
    f.close();
    GC.removeRange(f);
    free(cast(void*)f);
}

extern(C) SamHeader bamfile_get_header(BamFile* f) {
	return f.header;
}
