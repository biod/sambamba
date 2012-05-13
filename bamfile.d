module bamfile;

import bgzfrange;
import chunkinputstream;
import rangetransformer;
import samheader;
import reference;
import alignment;

import std.stream;
import std.system;
import std.algorithm : map;
import std.range : zip;
import std.zlib : uncompress, crc32, ZlibException;
import std.conv : to;
import std.exception : enforce;
import std.parallelism;

ubyte[] decompress(const BgzfBlock block) {
    auto uncompressed = uncompress(cast(void[])block.compressed_data, 
                                   cast(uint)block.input_size, -15);

    assert(block.input_size == uncompressed.length);
    assert(block.crc32 == crc32(0, uncompressed));

    return cast(ubyte[])uncompressed;
}


/**
  Represents BAM file
 */
struct BamFile {

    /**
      Constructor taking filename of BAM file to open,
      and optionally, task pool to use.
      
      Currently, opens the file read-only since library
      has no support for writing yet.
     */
    this(string filename, TaskPool task_pool = taskPool) {

        _filename = filename;
        _task_pool = task_pool;
        initializeStreams();

        auto magic = _bam.readString(4);
        
        enforce(magic == "BAM\1");

        readSamHeader();
        readReferenceSequencesInfo();

        // right after constructing, we are at the beginning
        //                           of the list of alignments
    }
    
    /*
       Get SAM header of file.
     */
    SamHeader header() @property {
        return _header;
    }

    /**
        Returns: information about reference sequences
     */
    ReferenceSequenceInfo[] reference_sequences() @property {
        return _reference_sequences;
    }

    /**
        Returns: range of alignments.
     */
    auto alignments() @property {
        if (_alignments_first_call) {
            // save rewind call
            _alignments_first_call = false;
        } else {
            rewind(); // if not the first call, need to rewind
        }
        return alignmentRange(_bam, _task_pool);
    }
    private bool _alignments_first_call = true;

    /**
        Seeks to the beginning of the list of alignments.
     */
    void rewind() {
        initializeStreams();
        _bam.readString(4); // skip magic
        int l_text;
        _bam.read(l_text);
        _bam.readString(l_text); // skip header
        int n_ref;
        _bam.read(n_ref);
        while (--n_ref > 0) {
            int l_name;
            _bam.read(l_name);
            _bam.readString(l_name);
            int l_ref;
            _bam.read(l_ref);
        } // skip reference sequences information
    }

    /*
       Closes underlying file stream
     */
    void close() {
        _file.close();
    }

private:
    string _filename;
    Stream _file;
    Stream _compressed_stream;
    BgzfRange _bgzf_range;
    Stream _bam;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    TaskPool _task_pool;

    // sets up the streams and ranges
    void initializeStreams() {
        
        _file = new BufferedFile(_filename);
        _compressed_stream = new EndianStream(_file, Endian.littleEndian);
        _bgzf_range = new BgzfRange(_compressed_stream);

        version(serial) {
            auto chunk_range = map!decompress(_bgzf_range);
        } else {
            /* TODO: tweak granularity */
            auto chunk_range = _task_pool.map!decompress(_bgzf_range, 25); 
        }
        
        auto decompressed_stream = makeChunkInputStream(chunk_range);
        _bam = new EndianStream(decompressed_stream, Endian.littleEndian); 
    }

    // initializes _header
    void readSamHeader() {
        int l_text;
        _bam.read(l_text);

        string text = to!string(_bam.readString(l_text));
        _header = SamHeader(text);
    }

    // initialize _reference_sequences
    void readReferenceSequencesInfo() {
        int n_ref;
        _bam.read(n_ref);
        _reference_sequences = new ReferenceSequenceInfo[n_ref];
        foreach (i; 0..n_ref) {
            _reference_sequences[i] = ReferenceSequenceInfo(_bam);
        }
    }
}
