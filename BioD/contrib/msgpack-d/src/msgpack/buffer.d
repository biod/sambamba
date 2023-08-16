// Written in the D programming language.

module msgpack.buffer;

//import std.traits;
import std.range;


version(Posix)
{
    import core.sys.posix.sys.uio : iovec;
}
else
{
    /**
     * from core.sys.posix.sys.uio.iovec for compatibility with posix.
     */
    struct iovec
    {
        void*  iov_base;
        size_t iov_len;
    }
}


/**
 * $(D RefBuffer) is a reference stored buffer for more efficient serialization
 *
 * Example:
 * -----
 * auto packer = packer(RefBuffer(16));  // threshold is 16
 *
 * // packs data
 *
 * writev(fd, cast(void*)packer.buffer.vector.ptr, packer.buffer.vector.length);
 * -----
 */
struct RefBuffer
{
  private:
    static struct Chunk
    {
        ubyte[] data;  // storing serialized value
        size_t  used;  // used size of data
    }

    immutable size_t Threshold;
    immutable size_t ChunkSize;

    // for putCopy
    Chunk[] chunks_;  // memory chunk for buffer
    size_t  index_;   // index for cunrrent chunk

    // for putRef
    iovec[] vecList_;  // reference to large data or copied data.


  public:
    /**
     * Constructs a buffer.
     *
     * Params:
     *  threshold = the threshold of writing value or stores reference.
     *  chunkSize = the default size of chunk for allocation.
     */
    @safe
    this(in size_t threshold, in size_t chunkSize = 8192)
    {
        Threshold = threshold;
        ChunkSize = chunkSize;

        chunks_.length = 1;
        chunks_[index_].data.length = chunkSize;
    }


    /**
     * Returns the buffer contents that excluding references.
     *
     * Returns:
     *  the non-contiguous copied contents.
     */
    @property @safe
    nothrow ubyte[] data()
    {
        ubyte[] result;

        foreach (ref chunk; chunks_)
            result ~= chunk.data[0..chunk.used];

        return result;
    }


    /**
     * Forwards to all buffer contents.
     *
     * Returns:
     *  the array of iovec struct that stores references.
     */
    @property @safe
    nothrow ref iovec[] vector() return
    {
        return vecList_;
    }


    /**
     * Writes the argument to buffer and stores the reference of writed content
     * if the argument size is smaller than threshold,
     * otherwise stores the reference of argument directly.
     *
     * Params:
     *  value = the content to write.
     */
    @safe
    void put(in ubyte value)
    {
        ubyte[1] values = [value];
        putCopy(values);
    }


    /// ditto
    @safe
    void put(in ubyte[] value)
    {
        if (value.length < Threshold)
            putCopy(value);
        else
            putRef(value);
    }


  private:
    /*
     * Stores the reference of $(D_PARAM value).
     *
     * Params:
     *  value = the content to write.
     */
    @trusted
    void putRef(in ubyte[] value)
    {
        vecList_.length += 1;
        vecList_[$ - 1]  = iovec(cast(void*)value.ptr, value.length);
    }


    /*
     * Writes $(D_PARAM value) to buffer and appends to its reference.
     *
     * Params:
     *  value = the contents to write.
     */
    @trusted
    void putCopy(const scope ubyte[] value)
    {
        /*
         * Helper for expanding new space.
         */
        void expand(in size_t size)
        {
            const newSize = size < ChunkSize ? ChunkSize : size;

            index_++;
            chunks_.length = 1;
            chunks_[index_].data.length = newSize;
        }

        const size = value.length;

        // lacks current chunk?
        if (chunks_[index_].data.length - chunks_[index_].used < size)
            expand(size);

        const base = chunks_[index_].used;                     // start index
        auto  data = chunks_[index_].data[base..base + size];  // chunk to write

        data[] = value[];
        chunks_[index_].used += size;

        // Optimization for avoiding iovec allocation.
        if (vecList_.length && data.ptr == (vecList_[$ - 1].iov_base +
                                            vecList_[$ - 1].iov_len))
            vecList_[$ - 1].iov_len += size;
        else
            putRef(data);
    }
}


unittest
{
    static assert(isOutputRange!(RefBuffer, ubyte) &&
                  isOutputRange!(RefBuffer, ubyte[]));

    auto buffer = RefBuffer(2, 4);

    ubyte[] tests = [1, 2];
    foreach (v; tests)
        buffer.put(v);
    buffer.put(tests);

    assert(buffer.data == tests, "putCopy failed");

    iovec[] vector = buffer.vector;
    ubyte[] result;

    assert(vector.length == 2, "Optimization failed");

    foreach (v; vector)
        result ~= (cast(ubyte*)v.iov_base)[0..v.iov_len];

    assert(result == tests ~ tests);
}
