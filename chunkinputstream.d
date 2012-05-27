module chunkinputstream;

import bgzfrange;
import virtualoffset;

import std.stream;
import std.range;

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
}

/**
  Class for turning range of DecompressedBgzfBlock objects
  into a proper InputStream
 */
final class ChunkInputStream(ChunkRange) : IChunkInputStream
{

    //static assert(ElementType!ChunkRange == DecompressedBgzfBlock);

    /// Construct class instance from range of chunks.
    this(ChunkRange range) 
    {
        _range = range;
        setupStream();
        readable = true;
        writeable = false;
        seekable = false;
    }

    override size_t readBlock(void* buffer, size_t size) {

        if (_range.empty()) {
            // no more chunks
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

private:
    ChunkRange _range;
    typeof(_range.front.start_offset) _start_offset;
    typeof(_range.front.decompressed_data) _buf;

    void setupStream() {

        if (_range.empty()) {
            readEOF = true;
            return;
        }

        auto tmp = _range.front; /// _range might be lazy, 
                                 /// so extra front() calls
                                 /// can cost a lot
        _start_offset = tmp.start_offset;
        _buf = tmp.decompressed_data;
        _len = _buf.length;

        if (_len == 0) {
            readEOF = true;
            return;
        }

        _cur = 0;

        return;
    }

    // We try to avoid ANY overhead in serial parts of the library,
    // and creating new MemoryStream for each block is expensive.
    // So needed parts are just copy-pasted from phobos/std/stream.d
    // with minor modifications.

    ulong _len;  // current data length
    ulong _cur;  // current position

    size_t readBlockHelper(void* buffer, size_t size) {
        ubyte* cbuf = cast(ubyte*) buffer;
        if (_len - _cur < size)
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
auto makeChunkInputStream(R)(R range) {
    return new ChunkInputStream!R(range);
}
