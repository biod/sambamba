import bgzfrange;
import baminputstream;
import rangetransformer;

import std.stdio;
import std.stream;
import std.system;
import std.algorithm : map;
import std.range : zip;
import std.zlib : uncompress, crc32, ZlibException;
import std.conv : to;

import utils.inputrangechunks;
import std.algorithm;
import std.array;
import std.functional;

ubyte[] decompress(const BgzfBlock block) {
    auto uncompressed = uncompress(cast(void[])block.compressed_data, 
                                   cast(uint)block.input_size, -15);

    assert(block.input_size == uncompressed.length);
    assert(block.crc32 == crc32(0, uncompressed));

    return cast(ubyte[])uncompressed;
}

auto chunkdecompress(BgzfBlock[] blocks) {
    return array(map!decompress(blocks));
}

int main(string[] args) {
    
    BufferedFile file;
    if (args.length == 2) {
        file = new BufferedFile(args[1]);
    } else {
        writeln(q{
usage: bamreader <filename>
                });
        return 0;
    }
    
    try {
        auto bgzf_range = chunkedInputRange(new BgzfRange(file), 150);
   
        import std.parallelism;
        auto pool = new TaskPool(2);
        scope(exit) pool.finish();

        auto chunk_range = joiner(rangeTransformer!chunkdecompress(bgzf_range, pool));

        auto stream = new EndianStream(
                          new BamInputStream!(typeof(chunk_range))(chunk_range),
                          Endian.littleEndian);

        auto magic = stream.readString(4);
        if (magic != "BAM\1") {
            writeln("not a BAM file!");
            writeln(magic);
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
