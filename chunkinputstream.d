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

        size_t read = _stream.readBlock(buffer, size);

        if (read < size) {
            _range.popFront(); // get next chunk
            setupStream();
        }

        if (read == 0) { // try with next chunk
            read = _stream.readBlock(buffer, size);
        }

        if (_stream.eof()) {
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
	MemoryStream _stream;

    private void setupStream() {

        if (!_range.empty()) {
            assert(_range.front.length > 0);
            _stream = new MemoryStream(_range.front);
            return;
        }
        
        readEOF = true;
    }
}

/**
   Returns: input stream wrapping given range of memory chunks
 */
auto makeChunkInputStream(R)(R range) {
    return new ChunkInputStream!R(range);
}
