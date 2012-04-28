import bgzfrange;
import baminputstream;

import std.stdio;
import std.stream;
import std.system;
import std.algorithm : map;
import std.range : zip;
import std.zlib : uncompress, crc32, ZlibException;
import std.conv : to;

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
        writeln("usage: bamreader <filename>\n       prints info about reference sequences: name, length, and\n       in how many alignment records its ID was referenced");
        return 0;
    }
    
    try {
        auto bgzf_range = new BgzfRange(file);
        auto chunk_range = map!decompress(bgzf_range);
        auto stream = new EndianStream(new BamInputStream!(typeof(chunk_range))(chunk_range),
                                       Endian.littleEndian);

        auto magic = stream.readString(4);
        if (magic != "BAM\1") {
            writeln("not a BAM file!");
            return 0;
        }

        int header_length;
        stream.read(header_length);

        stream.readString(header_length);

        int n_ref;
        stream.read(n_ref);

        struct ReferenceSequence {
            string name;
            int length;
        }

        ReferenceSequence[] refseqs;

        foreach (k; 0..n_ref) {
            int l_name;
            int l_ref;
            stream.read(l_name);
            string name = stream.readString(l_name).idup;
            stream.read(l_ref);
            
            refseqs ~= ReferenceSequence(name, l_ref);
        }
      
        uint[] occurrences = new uint[n_ref];
        
        while (!stream.eof()) {
            int block_size;
            stream.read(block_size);

            int ref_id;
            stream.read(ref_id);
            if (ref_id != -1 && ref_id < n_ref) {
                ++occurrences[ref_id];
            } 

            stream.readString(block_size - ref_id.sizeof);
        }

        foreach (seq; zip(refseqs, occurrences)) {
            writefln("%s (length %d): referenced %d times", 
                     seq[0].name, seq[0].length, seq[1]);
        }

    } catch (Exception e) {
        writeln(e);
        return 1;
    }
    return 0;
}
