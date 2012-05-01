module baminputstream;

import std.stream;

class BamInputStream(ChunkRange) : Stream {
    ChunkRange range;
	MemoryStream stream;

    private void setupStream() {

        if (!range.empty()) {
            assert(range.front.length > 0);
            stream = new MemoryStream(range.front);
            return;
        }
        
        readEOF = true;
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
            // no more chunks
            readEOF = true;
            return 0;
        }

        size_t read = stream.readBlock(buffer, size);

        if (read < size) {
            range.popFront(); // get next chunk
            setupStream();
        }

        if (read == 0) { // try with next chunk
            read = stream.readBlock(buffer, size);
        }

        if (stream.eof()) {
            range.popFront();

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
}
