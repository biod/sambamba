module chunkinputstream;

import bgzfrange;
import virtualoffset;

import std.stream;
import std.range;
import std.array;

/// Interface on top of Stream, which glues decompressed blocks
/// together and provides virtualTell() method for determining
/// virtual offset in the stream.
///
/// NOTE: the class has no virtualSeek method. The reason is
///       that an implementing class should have no direct access 
///       to any of underlying streams, (which provide chunks 
///       and their start virtual offsets), the reason being
///       the single responsibility principle. 
///       
///       If you need to do "virtualSeek", seek to coffset
///       in stream providing compressed blocks, create new
///       stream for decompressing, and skip uoffset bytes
///       in the decompressed stream. In general, you probably
///       will want to use different  decompression function
///       for each use case. (For serial iteration - without any 
///       caching, for random access - with caching; and maybe 
///       you'll want several types of caching.)
///
abstract class IChunkInputStream : Stream {
    /// Returns: current virtual offset
    ///
    /// (For details about what virtual offset is, 
    /// see bgzfrange module documentation.)
    ///
    abstract VirtualOffset virtualTell();

    /// Read a slice from chunk stream.
    abstract ubyte[] readSlice(size_t n);

    /// Average compression so far
    abstract float average_compression_ratio() @property const;

    /// Size of underlying BAM file (0 if not available).
    abstract size_t compressed_file_size() @property const;
}

/**
  Class for turning range of DecompressedBgzfBlock objects
  into a proper InputStream
 */
final class ChunkInputStream(ChunkRange) : IChunkInputStream
{

    static assert(is(ElementType!ChunkRange == DecompressedBgzfBlock));

    /// Construct class instance from range of chunks.
    this(ChunkRange range, size_t compressed_file_size=0) 
    {
        _range = range;
        setupStream();
        readable = true;
        writeable = false;
        seekable = false;

        _compressed_file_size = compressed_file_size;
    }

    /// Read a slice from chunk stream.
    ///
    /// Behaviour: when current chunk contains >= n bytes,
    ///            a chunk slice is returned. 
    ///            Otherwise, memory copying occurs.
    ubyte[] readSlice(size_t n) {
        if (_range.empty()) {
            _start_offset = _end_offset;
            _cur = 0;
            readEOF = true;
            return null;
        }

        if (_cur + n < _len) {
            // can return a slice, 
            // remaining in the current chunk
            auto slice = _buf[_cur .. _cur + n];
            _cur += n;
            return slice;
        } else if (_cur + n == _len) {
            // end of current chunk
            auto slice = _buf[_cur .. _len];
            _range.popFront();
            setupStream();
            return slice;
        } else {
            // this branch will be executed rarely
            // (wish the compiler could understand that...)
            auto slice = uninitializedArray!(ubyte[])(n);
            readExact(slice.ptr, n);
            return slice;
        }
    }

    override size_t readBlock(void* buffer, size_t size) {

        if (_range.empty()) {
            // no more chunks
            _start_offset = _end_offset;
            _cur = 0;
            readEOF = true;
            return 0;
        }

        size_t read = readBlockHelper(buffer, size);

        if (read < size) {
            _range.popFront(); // get next chunk
            setupStream();
        }

        if (read == 0) { // try with next chunk
            read = readBlockHelper(buffer, size);
        }

        if (_cur == _len) { // end of current chunk
            _range.popFront();

            // that will update readEOF flag if necessary
            setupStream(); 
        }

        return read;
    }

    override size_t writeBlock(const void* buffer, size_t size) {
        throw new WriteException("Stream is not writeable");
    }

    override ulong seek(long offset, SeekPos whence) {
        throw new SeekException("Stream is not seekable");
    }

    /// Returns: current virtual offset
    VirtualOffset virtualTell() {
		assert(_cur < (1<<16));
        return VirtualOffset(_start_offset, cast(ushort)_cur);
    }

    float average_compression_ratio() @property const {
        if (_total_compressed == 0) return 0.0;
        return cast(float)_total_uncompressed/_total_compressed;
    }

    size_t compressed_file_size() @property const {
        return _compressed_file_size;
    }

private:
    ChunkRange _range;
    typeof(_range.front.start_offset) _start_offset;
    typeof(_range.front.end_offset) _end_offset;
    typeof(_range.front.decompressed_data) _buf;

    size_t _compressed_file_size;

    size_t _len;  // current data length
    size_t _cur;  // current position

    size_t _total_compressed; // compressed bytes read so far
    size_t _total_uncompressed; // uncompressed size of blocks read so far

    void setupStream() {

        _cur = 0;

        if (_range.empty()) {
            _start_offset = _end_offset;
            readEOF = true;
            return;
        }

        auto tmp = _range.front; /// _range might be lazy, 
                                 /// so extra front() calls
                                 /// can cost a lot
        _start_offset = tmp.start_offset;
        _end_offset = tmp.end_offset;

        _total_compressed += _end_offset - _start_offset;

        _buf = tmp.decompressed_data;
        _len = _buf.length;

        _total_uncompressed += _len;

        if (_len == 0) {
            readEOF = true;
        }

        return;
    }

    size_t readBlockHelper(void* buffer, size_t size) {
        ubyte* cbuf = cast(ubyte*) buffer;
        if (size + _cur > _len)
          size = cast(size_t)(_len - _cur);
        ubyte[] ubuf = cast(ubyte[])_buf[cast(size_t)_cur .. cast(size_t)(_cur + size)];
        cbuf[0 .. size] = ubuf[];
        _cur += size;
        return size;
    }

}

/**
   Returns: input stream wrapping given range of memory chunks
 */
auto makeChunkInputStream(R)(R range, size_t compressed_file_size=0) {
    return new ChunkInputStream!R(range, compressed_file_size);
}
