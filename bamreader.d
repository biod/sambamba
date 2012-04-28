import bgzfrange;
import baminputstream;

import std.stdio : writeln;
import std.stream : BufferedFile;
import std.algorithm : map;
import std.zlib : uncompress, crc32, ZlibException;

ubyte[] decompress(BgzfBlock block) {
    auto uncompressed = uncompress(block.compressed_data, 
                                   block.input_size, -15);

    assert(block.input_size == uncompressed.length);
    assert(block.crc32 == crc32(0, uncompressed));

    return cast(ubyte[])uncompressed;
}

int main(string[] args) {
    BufferedFile file;

    if (args.length == 2) {
        file = new BufferedFile(args[1]);
    } else {
        writeln("usage: bamreader <filename>\n       prints SAM header");
        return 0;
    }
    
    try {
        auto bgzf_range = new BgzfRange(file);
        auto chunk_range = map!decompress(bgzf_range);
        auto stream = new BamInputStream!(typeof(chunk_range))(chunk_range);

        auto magic = stream.readString(4);
        if (magic != "BAM\1") {
            writeln("not a BAM file!");
            return 0;
        }

        int header_length;
        stream.read(header_length);

        writeln(stream.readString(header_length));
        
    } catch (ZlibException e) {
        writeln(e);
        return 1;
    }
    return 0;
}
