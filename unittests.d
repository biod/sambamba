/**
  Integration tests
 */

import samheader;
import bamfile;

unittest {
    import std.path;
    import std.stdio;

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
}

void main() {
}
