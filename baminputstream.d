module baminputstream;

import std.stream;

class BamInputStream(ChunkRange) : Stream {
    ChunkRange range;
    MemoryStream stream;

    private void setupStream() {
        if (!range.empty()) {
            stream = new MemoryStream(cast(ubyte[])range.front);
        }
    }

    public this(ChunkRange range) {
        this.range = range;
        setupStream();
        readable = true;
        writeable = false;
        seekable = false;
    }

    override size_t readBlock(void* buffer, size_t size) {
        if (range.empty()) {
            return 0; // EOF
        }
        size_t read = stream.readBlock(buffer, size);
        if (read < size) {
            range.popFront(); // get next chunk
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
}
