module chunkinputstream;

import std.stream;

/**
  Class for turning range of chunks (void[]/ubyte[] arrays)
  into a proper InputStream
 */
class ChunkInputStream(ChunkRange) : Stream {

    this(ChunkRange range) {
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

private:
    ChunkRange _range;
    typeof(_range.front) _buf;

    void setupStream() {

        if (_range.empty()) {
            readEOF = true;
            return;
        }

        _buf = _range.front;
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
