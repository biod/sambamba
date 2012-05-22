/**
  Integration tests
 */

import samheader;
import bamfile;
import bgzfrange;

import std.path;
import std.stdio;
import std.algorithm;
import std.conv;
import std.exception;

unittest {

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

    writeln("Testing alignment parsing...");
    fn = buildPath(dirName(__FILE__), "test", "data", "ex1_header.bam");
    bf = BamFile(fn);
    auto alignments = bf.alignments;
    auto read = alignments.front;
    assert(read.sequence == "CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA");
    assert(equal(map!"cast(char)(a + 33)"(read.phred_base_quality),
                "<<<7<<<;<<<<<<<<8;;<7;4<;<;;;;;94<;"));
    assert(bf.reference_sequences[read.ref_id].name == "chr1");
    assert(read.read_name == "EAS56_57:6:190:289:82");
    assert(read.flag == 69);
    assert(read.position == 99);
    assert(read.mapping_quality == 0);
    alignments.popFront();
    alignments.popFront();
    assert(alignments.front.cigar_string == "35M");

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
        read_name = read_name[1..$];
        string value = alignment.tags[tag.idup].to_sam;
        if (read_name != value) {
            writeln("tag: ", tag, "\tread_name: ", read_name, "\tvalue: ", value);
            writeln("value bam_typeid: ", alignment.tags[tag.idup].bam_typeid);
        }
        assert(read_name == value);
    }

    writeln("Testing exception handling...");
    fn = buildPath(dirName(__FILE__), "test", "data", "duplicated_block_size.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "no_block_size.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "wrong_extra_gzip_length.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "wrong_bc_subfield_length.bam");
    assertThrown!BgzfException(BamFile(fn));
    fn = buildPath(dirName(__FILE__), "test", "data", "corrupted_zlib_archive.bam");
    assertThrown!ZlibException(BamFile(fn));
}

void main() {
}
