module bamfile;

public import samheader;
public import reference;
public import alignment;
public import virtualoffset;
public import tagvalue;
public import alignmentrange;
import bgzfrange;
import chunkinputstream;
import randomaccessmanager;
import bai.read;
import utils.range;

import std.stream;
import std.system;
import std.stdio;
import std.algorithm : map, min;
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
    auto alignments(alias IteratePolicy=withoutOffsets)() @property {
		auto _decompressed_stream = getDecompressedAlignmentStream();
		return alignmentRange!IteratePolicy(_decompressed_stream);
    }

    /**
        Returns: range of all alignments in the file, calling $(D progressBarFunc)
                 for each alignment. 

        $(D progressBarFunc) will be called
        each time next alignment is read, with the argument being a number from [0.0, 1.0],
        which is estimated progress percentage.
    */
    auto alignmentsWithProgress(alias IteratePolicy=withoutOffsets)
        (void delegate(lazy float p) progressBarFunc) 
    {
		auto _decompressed_stream = getDecompressedAlignmentStream();
        auto alignments_with_offsets = alignmentRange!withOffsets(_decompressed_stream);

        static struct Result(alias IteratePolicy, R, S) {
            this(R range, S stream, void delegate(lazy float p) progressBarFunc) {
                _range = range;
                _stream = stream;
                _progress_bar_func = progressBarFunc;
            }

            static if (__traits(identifier, IteratePolicy) == "withOffsets") {
                auto front() @property {
                    return _range.front;
                } 
            } else static if (__traits(identifier, IteratePolicy) == "withoutOffsets") {
                auto front() @property {
                    return _range.front.alignment;
                }
            } else static assert(0, __traits(identifier, IteratePolicy));

            bool empty() @property {
                return _range.empty;
            }

            void popFront() {
                _bytes_read += _range.front.alignment.size_in_bytes;
                _range.popFront();
                if (_progress_bar_func !is null) {
                    _progress_bar_func(min(1.0, 
                        cast(float)_bytes_read / (_stream.compressed_file_size * 
                                                  _stream.average_compression_ratio)));
                }
            }

            private R _range;
            private S _stream;
            private size_t _bytes_read;
            private void delegate(lazy float p) _progress_bar_func;
        }

        return Result!(IteratePolicy, 
                       typeof(alignments_with_offsets),
                       typeof(_decompressed_stream))(alignments_with_offsets, 
                                                     _decompressed_stream,
                                                     progressBarFunc);
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

    RandomAccessManager _random_access_manager;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;
    int[string] _reference_sequence_dict; /// name -> index mapping

    TaskPool _task_pool;
    size_t buffer_size = 1048576; // buffer size to be used for I/O

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
   
         if (compressed_stream.seekable) {
            return makeChunkInputStream(chunk_range, cast(size_t)file.size);
        } else {
            return makeChunkInputStream(chunk_range);
        }
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
    
        auto sz = compressed_stream.seekable ? compressed_stream.size : 0;
		auto stream = makeChunkInputStream(chunk_range, cast(size_t)sz);
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
