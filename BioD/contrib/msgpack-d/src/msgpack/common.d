// Written in the D programming language.

module msgpack.common;

import msgpack.attribute;

import std.typetuple; // will use std.meta
import std.traits;


// for Converting Endian using ntohs and ntohl;
version(Windows)
{
    package import core.sys.windows.winsock2;
}
else
{
    package import core.sys.posix.arpa.inet;
}


version(EnableReal)
{
    package enum EnableReal = true;
}
else
{
    package enum EnableReal = false;
}


static if (real.sizeof == double.sizeof) {
    // for 80bit real inter-operation on non-x86 CPU
    version = NonX86;

    package import std.numeric;
}


@trusted:
public:


/**
 * $(D ExtValue) is a $(D MessagePack) Extended value representation.
 * The application is responsible for correctly interpreting $(D data) according
 *  to the type described by $(D type).
 */
struct ExtValue
{
    byte type;    /// An integer 0-127 with application-defined meaning
    ubyte[] data; /// The raw bytes
}


/**
 * MessagePack type-information format
 *
 * See_Also:
 *  $(LINK2 http://redmine.msgpack.org/projects/msgpack/wiki/FormatSpec, MessagePack Specificaton)
 */
enum Format : ubyte
{
    // unsinged integer
    UINT8  = 0xcc,  // ubyte
    UINT16 = 0xcd,  // ushort
    UINT32 = 0xce,  // uint
    UINT64 = 0xcf,  // ulong

    // signed integer
    INT8  = 0xd0,   // byte
    INT16 = 0xd1,   // short
    INT32 = 0xd2,   // int
    INT64 = 0xd3,   // long

    // floating point
    FLOAT  = 0xca,  // float
    DOUBLE = 0xcb,  // double

    // raw byte
    RAW   = 0xa0,
    RAW16 = 0xda,
    RAW32 = 0xdb,

    // bin type
    BIN8  = 0xc4,
    BIN16 = 0xc5,
    BIN32 = 0xc6,

    // ext type
    EXT   = 0xd4,  // fixext 1/2/4/8/16
    EXT8  = 0xc7,
    EXT16 = 0xc8,
    EXT32 = 0xc9,

    // str type
    STR8  = 0xd9,
    //STR16 = 0xda,
    //STR32 = 0xdb,

    // array
    ARRAY   = 0x90,
    ARRAY16 = 0xdc,
    ARRAY32 = 0xdd,

    // map
    MAP   = 0x80,
    MAP16 = 0xde,
    MAP32 = 0xdf,

    // other
    NIL   = 0xc0,   // null
    TRUE  = 0xc3,
    FALSE = 0xc2,

    // real (This format is D only!)
    REAL = 0xd4
}


package:


/**
 * For float type serialization / deserialization
 */
union _f
{
    float f;
    uint  i;
}


/**
 * For double type serialization / deserialization
 */
union _d
{
    double f;
    ulong  i;
}


/**
 * For real type serialization / deserialization
 *
 * 80-bit real is padded to 12 bytes(Linux) and 16 bytes(Mac).
 * http://lists.puremagic.com/pipermail/digitalmars-d/2010-June/077394.html
 */
union _r
{
    real f;

    struct
    {
        ulong  fraction;
        ushort exponent;  // includes sign
    }
}

enum RealSize = 10;  // Real size is 80bit


/**
 * Detects whether $(D_PARAM T) is a built-in byte type.
 */
template isByte(T)
{
    enum isByte = staticIndexOf!(Unqual!T, byte, ubyte) >= 0;
}


unittest
{
    static assert(isByte!(byte));
    static assert(isByte!(const(byte)));
    static assert(isByte!(ubyte));
    static assert(isByte!(immutable(ubyte)));
    static assert(!isByte!(short));
    static assert(!isByte!(char));
    static assert(!isByte!(string));
}


/**
 * Gets asterisk string from pointer type
 */
template AsteriskOf(T)
{
    static if (is(T P == U*, U))
        enum AsteriskOf = "*" ~ AsteriskOf!U;
    else
        enum AsteriskOf = "";
}


/**
 * Get the number of member to serialize.
 */
template SerializingMemberNumbers(Classes...)
{
    static if (Classes.length == 0)
        enum SerializingMemberNumbers = 0;
    else
        enum SerializingMemberNumbers = Filter!(isPackedField, Classes[0].tupleof).length + SerializingMemberNumbers!(Classes[1..$]);
}


/**
 * Get derived classes with serialization-order
 */
template SerializingClasses(T)
{
    // There is no information in Object type. Currently disable Object serialization.
    static if (is(T == Object))
        static assert(false, "Object type serialization doesn't support yet. Please define toMsgpack/fromMsgpack and use cast");
    else
        alias TypeTuple!(Reverse!(Erase!(Object, BaseClassesTuple!(T))), T) SerializingClasses;
}


/**
 * Get a field name of class or struct.
 */
template getFieldName(Type, size_t i)
{
    import std.conv : text;

    static assert((is(Unqual!Type == class) || is(Unqual!Type == struct)), "Type must be class or struct: type = " ~ Type.stringof);
    static assert(i < Type.tupleof.length, text(Type.stringof, " has ", Type.tupleof.length, " attributes: given index = ", i));

    enum getFieldName = __traits(identifier, Type.tupleof[i]);
}


version (LittleEndian)
{
    /*
     * Converts $(value) to different Endian.
     *
     * Params:
     *  value = the LittleEndian value to convert.
     *
     * Returns:
     *  the converted value.
     */
    @trusted
    ushort convertEndianTo(size_t Bit, T)(in T value) if (Bit == 16)
    {
        return ntohs(cast(ushort)value);
    }


    // ditto
    @trusted
    uint convertEndianTo(size_t Bit, T)(in T value) if (Bit == 32)
    {
        return ntohl(cast(uint)value);
    }


    // ditto
    @trusted
    ulong convertEndianTo(size_t Bit, T)(in T value) if (Bit == 64)
    {
        // dmd has convert function?
        return ((((cast(ulong)value) << 56) & 0xff00000000000000UL) |
                (((cast(ulong)value) << 40) & 0x00ff000000000000UL) |
                (((cast(ulong)value) << 24) & 0x0000ff0000000000UL) |
                (((cast(ulong)value) <<  8) & 0x000000ff00000000UL) |
                (((cast(ulong)value) >>  8) & 0x00000000ff000000UL) |
                (((cast(ulong)value) >> 24) & 0x0000000000ff0000UL) |
                (((cast(ulong)value) >> 40) & 0x000000000000ff00UL) |
                (((cast(ulong)value) >> 56) & 0x00000000000000ffUL));
    }


    unittest
    {
        assert(convertEndianTo!16(0x0123)             == 0x2301);
        assert(convertEndianTo!32(0x01234567)         == 0x67452301);
        assert(convertEndianTo!64(0x0123456789abcdef) == 0xefcdab8967452301);
    }


    /*
     * Comapatible for BigEndian environment.
     */
    ubyte take8from(size_t bit = 8, T)(T value)
    {
        static if (bit == 8 || bit == 16 || bit == 32 || bit == 64)
            return (cast(ubyte*)&value)[0];
        else
            static assert(false, bit.stringof ~ " is not support bit width.");
    }


    unittest
    {
        foreach (Integer; TypeTuple!(ubyte, ushort, uint, ulong)) {
            assert(take8from!8 (cast(Integer)0x01)               == 0x01);
            assert(take8from!16(cast(Integer)0x0123)             == 0x23);
            assert(take8from!32(cast(Integer)0x01234567)         == 0x67);
            assert(take8from!64(cast(Integer)0x0123456789abcdef) == 0xef);
        }
    }
}
else
{
    /*
     * Comapatible for LittleEndian environment.
     */
    @safe
    ushort convertEndianTo(size_t Bit, T)(in T value) if (Bit == 16)
    {
        return cast(ushort)value;
    }


    // ditto
    @safe
    uint convertEndianTo(size_t Bit, T)(in T value) if (Bit == 32)
    {
        return cast(uint)value;
    }


    // ditto
    @safe
    ulong convertEndianTo(size_t Bit, T)(in T value) if (Bit == 64)
    {
        return cast(ulong)value;
    }


    unittest
    {
        assert(convertEndianTo!16(0x0123)       == 0x0123);
        assert(convertEndianTo!32(0x01234567)   == 0x01234567);
        assert(convertEndianTo!64(0x0123456789) == 0x0123456789);
    }


    /*
     * Takes 8bit from $(D_PARAM value)
     *
     * Params:
     *  value = the content to take.
     *
     * Returns:
     *  the 8bit value corresponding $(D_PARAM bit) width.
     */
    ubyte take8from(size_t bit = 8, T)(T value)
    {
        static if (bit == 8)
            return (cast(ubyte*)&value)[0];
        else static if (bit == 16)
            return (cast(ubyte*)&value)[1];
        else static if (bit == 32)
            return (cast(ubyte*)&value)[3];
        else static if (bit == 64)
            return (cast(ubyte*)&value)[7];
        else
            static assert(false, bit.stringof ~ " is not support bit width.");
    }


    unittest
    {
        foreach (Integer; TypeTuple!(ubyte, ushort, uint, ulong)) {
            assert(take8from!8 (cast(Integer)0x01)               == 0x01);
            assert(take8from!16(cast(Integer)0x0123)             == 0x23);
            assert(take8from!32(cast(Integer)0x01234567)         == 0x67);
            assert(take8from!64(cast(Integer)0x0123456789abcdef) == 0xef);
        }
    }
}


/*
 * Loads $(D_PARAM T) type value from $(D_PARAM buffer).
 *
 * Params:
 *  buffer = the serialized contents.
 *
 * Returns:
 *  the Endian-converted value.
 */
T load16To(T)(ubyte[] buffer)
{
    return cast(T)(convertEndianTo!16(*cast(ushort*)buffer.ptr));
}


// ditto
T load32To(T)(ubyte[] buffer)
{
    return cast(T)(convertEndianTo!32(*cast(uint*)buffer.ptr));
}


// ditto
T load64To(T)(ubyte[] buffer)
{
    return cast(T)(convertEndianTo!64(*cast(ulong*)buffer.ptr));
}


version (D_Ddoc)
{
    /**
     * Internal buffer and related operations for Unpacker
     *
     * Following Unpackers mixin this template. So, Unpacker can use following methods.
     *
     * -----
     * //buffer image:
     * +-------------------------------------------+
     * | [object] | [obj | unparsed... | unused... |
     * +-------------------------------------------+
     *            ^ offset
     *                   ^ current
     *                                 ^ used
     *                                             ^ buffer.length
     * -----
     *
     * This mixin template is a private.
     */
    mixin template InternalBuffer()
    {
      private:
        ubyte[] buffer_;  // internal buffer
        size_t  used_;    // index that buffer cosumed
        size_t  offset_;  // index that buffer parsed
        size_t  parsed_;  // total size of parsed message
        bool    hasRaw_;  // indicates whether Raw object has been deserialized


      public:
        /**
         * Forwards to internal buffer.
         *
         * Returns:
         *  the reference of internal buffer.
         */
        @property @safe
        nothrow ubyte[] buffer();


        /**
         * Fills internal buffer with $(D_PARAM target).
         *
         * Params:
         *  target = new serialized buffer to deserialize.
         */
        @safe void feed(in ubyte[] target);


        /**
         * Consumes buffer. This method is helper for buffer property.
         * You must use this method if you write bytes to buffer directly.
         *
         * Params:
         *  size = the number of consuming.
         */
        @safe
        nothrow void bufferConsumed(in size_t size);


        /**
         * Removes unparsed buffer.
         */
        @safe
        nothrow void removeUnparsed();


        /**
         * Returns:
         *  the total size including unparsed buffer size.
         */
        @property @safe
        nothrow size_t size() const;


        /**
         * Returns:
         *  the parsed size of buffer.
         */
        @property @safe
        nothrow size_t parsedSize() const;


        /**
         * Returns:
         *  the unparsed size of buffer.
         */
        @property @safe
        nothrow size_t unparsedSize() const;


    private:
        @safe
        void initializeBuffer(in ubyte[] target, in size_t bufferSize = 8192);
    }
}
else
{
    mixin template InternalBuffer()
    {
      private:
        ubyte[] buffer_;  // internal buffer
        size_t  used_;    // index that buffer cosumed
        size_t  offset_;  // index that buffer parsed
        size_t  parsed_;  // total size of parsed message
        bool    hasRaw_;  // indicates whether Raw object has been deserialized


      public:
        @property @safe
        nothrow ubyte[] buffer()
        {
            return buffer_;
        }


        @safe
        void feed(in ubyte[] target)
        in
        {
            assert(target.length);
        }
        do
        {
            /*
             * Expands internal buffer.
             *
             * Params:
             *  size = new buffer size to append.
             */
            void expandBuffer(in size_t size)
            {
                // rewinds buffer(completed deserialization)
                if (used_ == offset_ && !hasRaw_) {
                    used_ =  offset_ = 0;

                    if (buffer_.length < size)
                        buffer_.length = size;

                    return;
                }

                // deserializing state is mid-flow(buffer has non-parsed data yet)
                auto unparsed = buffer_[offset_..used_];
                auto restSize = buffer_.length - used_ + offset_;
                auto newSize  = size > restSize ? unparsedSize + size : buffer_.length;

                if (hasRaw_) {
                    hasRaw_ = false;
                    buffer_ = new ubyte[](newSize);
                } else {
                    buffer_.length = newSize;

                    // avoids overlapping copy
                    auto area = buffer_[0..unparsedSize];
                    unparsed  = area.overlap(unparsed) ? unparsed.dup : unparsed;
                }

                buffer_[0..unparsedSize] = unparsed[];
                used_   = unparsedSize;
                offset_ = 0;
            }

            const size = target.length;

            // lacks current buffer?
            if (buffer_.length - used_ < size)
                expandBuffer(size);

            buffer_[used_..used_ + size] = target[];
            used_ += size;
        }


        @safe
        nothrow void bufferConsumed(in size_t size)
        {
            if (used_ + size > buffer_.length)
                used_ = buffer_.length;
            else
                used_ += size;
        }


        @safe
        nothrow void removeUnparsed()
        {
            used_ = offset_;
        }


        @property @safe
        nothrow size_t size() const
        {
            return parsed_ - offset_ + used_;
        }


        @property @safe
        nothrow size_t parsedSize() const
        {
            return parsed_;
        }


        @property @safe
        nothrow size_t unparsedSize() const
        {
            return used_ - offset_;
        }


      private:
        @safe
        nothrow void initializeBuffer(in ubyte[] target, in size_t bufferSize = 8192)
        {
            const size = target.length;

            buffer_ = new ubyte[](size > bufferSize ? size : bufferSize);
            used_   = size;
            buffer_[0..size] = target[];
        }
    }
}
