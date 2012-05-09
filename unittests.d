/**
  Integration tests
 */

import samheader;
import bamfile;

unittest {
    import std.path;
    import std.stdio;
	import std.algorithm;
	import std.conv;

    writeln("Testing extracting SAM header...");
    auto fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    auto bf = BamFile(fn);
    assert(bf.header.format_version == "1.3");
    assert(bf.header.sorting_order == "coordinate");
    assert(bf.header.sq_lines.length == 2);
    assert(bf.header.sq_lines[0].sequence_name == "chr1");
    assert(bf.header.sq_lines[1].sequence_length == 1584);

    fn = buildPath(dirName(__FILE__), "test", "data", "bins.bam");
    bf = BamFile(fn);
    assert(bf.header.sorting_order == "unknown");
    assert(bf.header.sq_lines.length == 3);
    assert(bf.header.rg_lines.length == 0);
    assert(bf.header.sq_lines[2].sequence_name == "large");
    assert(bf.header.sq_lines[1].sequence_length == 65536);

	writeln("Testing tag parsing...");
	fn = buildPath(dirName(__FILE__), "test", "data", "tags.bam");
	bf = BamFile(fn);
	foreach (alignment; bf.alignments) {
		auto read_name = alignment.read_name;
		assert(read_name[0..4] == "tag_");
		char[] tag;
		read_name = read_name[4..$];
		while (read_name[0] != ':') {
			tag ~= read_name[0];
			read_name = read_name[1..$];
		}
		string expected;
		bool array = false;
		if (read_name[1] == 'B') {
			expected = read_name[5..$];
			array = true;
		} else {
			expected = read_name[3..$];
		}
		assert(alignment.tags.keys.length == 1);
		string value = to!string(alignment.tags[tag.idup]);
		if (array) {
			value = value[1..$-1]; // strip []
			value = to!string(filter!"a != ' '"(value));
		}
		if (expected != value) {
			writeln(tag, read_name, expected, value);
		}
		assert(expected == value);
	}
}

void main() {
}
