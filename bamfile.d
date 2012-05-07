module bamfile;

import bgzfrange;
import chunkinputstream;
import rangetransformer;
import samheader;

import std.stream;
import std.system;
import std.algorithm : map;
import std.range : zip;
import std.zlib : uncompress, crc32, ZlibException;
import std.conv : to;
import std.exception : enforce;

/**
  Represents BAM file
 */
struct BamFile {

    /**
      Constructor taking filename of BAM file to open.
      
      Currently, opens the file read-only since library
      has no support for writing yet.
     */
    this(string filename) {

        ubyte[] decompress(const BgzfBlock block) {
            auto uncompressed = uncompress(cast(void[])block.compressed_data, 
                                           cast(uint)block.input_size, -15);

            assert(block.input_size == uncompressed.length);
            assert(block.crc32 == crc32(0, uncompressed));

            return cast(ubyte[])uncompressed;
        }

        _file = new BufferedFile(filename);
        _compressed_stream = new EndianStream(_file, Endian.littleEndian);
        _bgzf_range = new BgzfRange(_compressed_stream);

        auto chunk_range = map!decompress(_bgzf_range); 
        
        auto decompressed_stream = makeChunkInputStream(chunk_range);
        _bam = new EndianStream(decompressed_stream, Endian.littleEndian); 

        auto magic = _bam.readString(4);
        
        enforce(magic == "BAM\1");

        readSamHeader();
    }
    
    /*
       Get SAM header of file.
     */
    SamHeader header() @property {
        return _header;
    }

    /*
       Close underlying file stream
     */
    void close() {
        _file.close();
    }

private:
    Stream _file;
    Stream _compressed_stream;
    BgzfRange _bgzf_range;
    Stream _bam;

    SamHeader _header;

    // initializes _header
    void readSamHeader() {
        int header_len;
        _bam.read(header_len);

        string header_contents = to!string(_bam.readString(header_len));
        _header = SamHeader(header_contents);
    }
}
