module bamfile;

public import samheader;
public import reference;
public import alignment;
public import virtualoffset;
public import tagvalue;
import alignmentrange;
import bgzfrange;
import chunkinputstream;
import randomaccessmanager;
import bai.read;
import utils.range;

import std.stream;
import std.system;
import std.stdio;
import std.algorithm : map;
import std.range : zip;
import std.conv : to;
import std.exception : enforce;
import std.parallelism;
import std.array : uninitializedArray;

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
      
        try {
            _bai_file = BaiFile(filename);
            _random_access_manager = new RandomAccessManager(_filename, _bai_file);
        } catch (Exception e) {
            // TODO: logging levels
            // stderr.writeln("Couldn't find index file: ", e.msg);
            _random_access_manager = new RandomAccessManager(_filename);
        }

        auto magic = _bam.readString(4);
        
        enforce(magic == "BAM\1", "Invalid file format: expected BAM\\1");

        readSamHeader();
        readReferenceSequencesInfo();

        // right after constructing, we are at the beginning
        //                           of the list of alignments

		_alignments_start_voffset = _decompressed_stream.virtualTell();
    }
  
    /// True if associated BAI file was found
    bool has_index() @property {
        return _random_access_manager.found_index_file;
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
        Returns: range of all alignments in the file.

        However, using several ranges is not recommended since it can hurt
		disk access performance.
     */
    auto alignments() @property {
		auto _decompressed_stream = getDecompressedAlignmentStream();
		return alignmentRange(_decompressed_stream);
    }

	/**
	  	Returns: range of all alignments in the file together with their
		start/end virtual offsets (alignmentrange.AlignmentBlock structs)
	 */
    auto alignmentsWithVirtualOffsets() {
		auto _decompressed_stream = getDecompressedAlignmentStream();
        return alignmentRangeWithOffsets(_decompressed_stream);
    }

    /**
      Get an alignment at a given virtual offset.
     */
    Alignment getAlignmentAt(VirtualOffset offset) {
        return _random_access_manager.getAlignmentAt(offset);
    }

    /**
      Returns reference sequence with id $(D ref_id).
     */
    auto reference(int ref_id) {
        enforce(ref_id < _reference_sequences.length, "Invalid reference index");
        return ReferenceSequence(_random_access_manager, 
                                 ref_id,
                                 _reference_sequences[ref_id]);
    }

    /**
      Returns reference sequence named $(D ref_name).
     */
    auto opIndex(string ref_name) {
        enforce(hasReference(ref_name), "Reference with name " ~ ref_name ~ " does not exist");
        auto ref_id = _reference_sequence_dict[ref_name];
        return reference(ref_id);
    }

    /**
      Check if reference named $(D ref_name) is presented in BAM header.
     */
    bool hasReference(string ref_name) {
        return null != (ref_name in _reference_sequence_dict);
    }

    /**
      Set buffer size for I/O operations.
     */
    void setBufferSize(size_t buffer_size) {
        this.buffer_size = buffer_size;
    }

private:
    
    string _filename;
    IChunkInputStream _decompressed_stream;
    Stream _bam;

	// Virtual offset at which alignment records start.
	VirtualOffset _alignments_start_voffset;

    BaiFile _bai_file; /// provides access to index file

    typeof(alignmentRange(_decompressed_stream)) _alignment_range;
    typeof(alignmentRangeWithOffsets(_decompressed_stream)) _alignment_range_with_offsets;

    RandomAccessManager _random_access_manager;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict; /// name -> index mapping

    TaskPool _task_pool;
    size_t buffer_size = 8192; // buffer size to be used for I/O

	// get decompressed stream out of compressed BAM file
	IChunkInputStream getDecompressedStream() {
		auto file = new BufferedFile(_filename, FileMode.In, buffer_size);
        auto compressed_stream = new EndianStream(file, Endian.littleEndian);
        auto bgzf_range = BgzfRange(compressed_stream);

        version(serial) {
            auto chunk_range = map!decompressBgzfBlock(bgzf_range);
        } else {
            auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 25);
        }
    	
		return makeChunkInputStream(chunk_range);
	}

	// get decompressed stream starting from the first alignment record
	IChunkInputStream getDecompressedAlignmentStream() {
		enforce(_alignments_start_voffset != 0UL);

		auto file = new BufferedFile(_filename, FileMode.In, buffer_size);
        auto compressed_stream = new EndianStream(file, Endian.littleEndian);
		compressed_stream.seekCur(_alignments_start_voffset.coffset);
        auto bgzf_range = BgzfRange(compressed_stream);

        version(serial) {
            auto chunk_range = map!decompressBgzfBlock(bgzf_range);
        } else {
            auto chunk_range = _task_pool.map!decompressBgzfBlock(bgzf_range, 25);
        }
    	
		auto stream = makeChunkInputStream(chunk_range);
		stream.readString(_alignments_start_voffset.uoffset);
		return stream;
	}

    // sets up the streams and ranges
    void initializeStreams() {
        
		_decompressed_stream = getDecompressedStream();
        _bam = new EndianStream(_decompressed_stream, Endian.littleEndian); 
    }

    // initializes _header
    void readSamHeader() {
        int l_text;
        _bam.read(l_text);

        string text = to!string(_bam.readString(l_text));
        _header = new SamHeader(text);
    }

    // initialize _reference_sequences
    void readReferenceSequencesInfo() {
        int n_ref;
        _bam.read(n_ref);
        _reference_sequences = new ReferenceSequenceInfo[n_ref];
        foreach (i; 0..n_ref) {
            _reference_sequences[i] = ReferenceSequenceInfo(_bam);

            // provide mapping Name -> Index
            _reference_sequence_dict[_reference_sequences[i].name] = i;
        }
    }
}
