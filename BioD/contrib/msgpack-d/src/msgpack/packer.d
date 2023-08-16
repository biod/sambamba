module msgpack.packer;

import msgpack.common;
import msgpack.attribute;
import msgpack.exception;

import std.array;
import std.exception;
import std.range;
import std.stdio;
import std.traits;
import std.typecons;
import std.typetuple;
import std.container;


/**
 * $(D Packer) is a $(D MessagePack) serializer
 *
 * Example:
 * -----
 * auto packer = packer(Appender!(ubyte[])());
 *
 * packer.packArray(false, 100, 1e-10, null);
 *
 * stdout.rawWrite(packer.stream.data);
 * -----
 *
 * NOTE:
 *  Current implementation can't deal with a circular reference.
 *  If you try to serialize a object that has circular reference, runtime raises 'Stack Overflow'.
 */
struct PackerImpl(Stream) if (isOutputRange!(Stream, ubyte) && isOutputRange!(Stream, ubyte[]))
{
  private:
    static @system
    {
        alias void delegate(ref PackerImpl, void*) PackHandler;
        PackHandler[TypeInfo] packHandlers;

        public void registerHandler(T, alias Handler)()
        {
            packHandlers[typeid(T)] = delegate(ref PackerImpl packer, void* obj) {
                Handler(packer, *cast(T*)obj);
            };
        }

        public void register(T)()
        {
            packHandlers[typeid(T)] = delegate(ref PackerImpl packer, void* obj) {
                packer.packObject(*cast(T*)obj);
            };
        }
    }

    enum size_t Offset = 1;  // type-information offset

    Stream                   stream_;  // the stream to write
    ubyte[Offset + RealSize] store_;   // stores serialized value
    bool                     withFieldName_;


  public:
    /**
     * Constructs a packer with $(D_PARAM stream).
     *
     * Params:
     *  stream        = the stream to write.
     *  withFieldName = serialize class / struct with field name
     */
    this(Stream stream, bool withFieldName = false)
    {
        stream_        = stream;
        withFieldName_ = withFieldName;
    }


    /**
     * Constructs a packer with $(D_PARAM withFieldName).
     *
     * Params:
     *  withFieldName = serialize class / struct with field name
     */
    this(bool withFieldName)
    {
        withFieldName_ = withFieldName;
    }


    /**
     * Forwards to stream.
     *
     * Returns:
     *  the stream.
     */
    @property @safe
    nothrow ref Stream stream()
    {
        return stream_;
    }


    /**
     * Serializes argument and writes to stream.
     *
     * If the argument is the pointer type, dereferences the pointer and serializes pointed value.
     * -----
     * int  a = 10;
     * int* b = &b;
     *
     * packer.pack(b);  // serializes 10, not address of a
     * -----
     * Serializes nil if the argument of nullable type is null.
     *
     * NOTE:
     *  MessagePack doesn't define $(D_KEYWORD real) type format.
     *  Don't serialize $(D_KEYWORD real) if you communicate with other languages.
     *  Transfer $(D_KEYWORD double) serialization if $(D_KEYWORD real) on your environment equals $(D_KEYWORD double).
     *
     * Params:
     *  value = the content to serialize.
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref PackerImpl pack(T)(in T value) if (is(Unqual!T == bool))
    {
        if (value)
            stream_.put(Format.TRUE);
        else
            stream_.put(Format.FALSE);

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (isUnsigned!T && !is(Unqual!T == enum))
    {
        // ulong < ulong is slower than uint < uint
        static if (!is(Unqual!T  == ulong)) {
            enum Bits = T.sizeof * 8;

            if (value < (1 << 8)) {
                if (value < (1 << 7)) {
                    // fixnum
                    stream_.put(take8from!Bits(value));
                } else {
                    // uint 8
                    store_[0] = Format.UINT8;
                    store_[1] = take8from!Bits(value);
                    stream_.put(store_[0..Offset + ubyte.sizeof]);
                }
            } else {
                if (value < (1 << 16)) {
                    // uint 16
                    const temp = convertEndianTo!16(value);

                    store_[0] = Format.UINT16;
                    *cast(ushort*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + ushort.sizeof]);
                } else {
                    // uint 32
                    const temp = convertEndianTo!32(value);

                    store_[0] = Format.UINT32;
                    *cast(uint*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + uint.sizeof]);
                }
            }
        } else {
            if (value < (1UL << 8)) {
                if (value < (1UL << 7)) {
                    // fixnum
                    stream_.put(take8from!64(value));
                } else {
                    // uint 8
                    store_[0] = Format.UINT8;
                    store_[1] = take8from!64(value);
                    stream_.put(store_[0..Offset + ubyte.sizeof]);
                }
            } else {
                if (value < (1UL << 16)) {
                    // uint 16
                    const temp = convertEndianTo!16(value);

                    store_[0] = Format.UINT16;
                    *cast(ushort*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + ushort.sizeof]);
                } else if (value < (1UL << 32)){
                    // uint 32
                    const temp = convertEndianTo!32(value);

                    store_[0] = Format.UINT32;
                    *cast(uint*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + uint.sizeof]);
                } else {
                    // uint 64
                    const temp = convertEndianTo!64(value);

                    store_[0] = Format.UINT64;
                    *cast(ulong*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + ulong.sizeof]);
                }
            }
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (isSigned!T && isIntegral!T && !is(Unqual!T == enum))
    {
        // long < long is slower than int < int
        static if (!is(Unqual!T == long)) {
            enum Bits = T.sizeof * 8;

            if (value < -(1 << 5)) {
                if (value < -(1 << 15)) {
                    // int 32
                    const temp = convertEndianTo!32(value);

                    store_[0] = Format.INT32;
                    *cast(int*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + int.sizeof]);
                } else if (value < -(1 << 7)) {
                    // int 16
                    const temp = convertEndianTo!16(value);

                    store_[0] = Format.INT16;
                    *cast(short*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + short.sizeof]);
                } else {
                    // int 8
                    store_[0] = Format.INT8;
                    store_[1] = take8from!Bits(value);
                    stream_.put(store_[0..Offset + byte.sizeof]);
                }
            } else if (value < (1 << 7)) {
                // fixnum
                stream_.put(take8from!Bits(value));
            } else {
                if (value < (1 << 8)) {
                    // uint 8
                    store_[0] = Format.UINT8;
                    store_[1] = take8from!Bits(value);
                    stream_.put(store_[0..Offset + ubyte.sizeof]);
                } else if (value < (1 << 16)) {
                    // uint 16
                    const temp = convertEndianTo!16(value);

                    store_[0] = Format.UINT16;
                    *cast(ushort*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + ushort.sizeof]);
                } else {
                    // uint 32
                    const temp = convertEndianTo!32(value);

                    store_[0] = Format.UINT32;
                    *cast(uint*)&store_[Offset] = temp;
                    stream_.put(store_[0..Offset + uint.sizeof]);
                }
            }
        } else {
            if (value < -(1L << 5)) {
                if (value < -(1L << 15)) {
                    if (value < -(1L << 31)) {
                        // int 64
                        const temp = convertEndianTo!64(value);

                        store_[0] = Format.INT64;
                        *cast(long*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + long.sizeof]);
                    } else {
                        // int 32
                        const temp = convertEndianTo!32(value);

                        store_[0] = Format.INT32;
                        *cast(int*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + int.sizeof]);
                    }
                } else {
                    if (value < -(1L << 7)) {
                        // int 16
                        const temp = convertEndianTo!16(value);

                        store_[0] = Format.INT16;
                        *cast(short*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + short.sizeof]);
                    } else {
                        // int 8
                        store_[0] = Format.INT8;
                        store_[1] = take8from!64(value);
                        stream_.put(store_[0..Offset + byte.sizeof]);
                    }
                }
            } else if (value < (1L << 7)) {
                // fixnum
                stream_.put(take8from!64(value));
            } else {
                if (value < (1L << 16)) {
                    if (value < (1L << 8)) {
                        // uint 8
                        store_[0] = Format.UINT8;
                        store_[1] = take8from!64(value);
                        stream_.put(store_[0..Offset + ubyte.sizeof]);
                    } else {
                        // uint 16
                        const temp = convertEndianTo!16(value);

                        store_[0] = Format.UINT16;
                        *cast(ushort*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + ushort.sizeof]);
                    }
                } else {
                    if (value < (1L << 32)) {
                        // uint 32
                        const temp = convertEndianTo!32(value);

                        store_[0] = Format.UINT32;
                        *cast(uint*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + uint.sizeof]);
                    } else {
                        // uint 64
                        const temp = convertEndianTo!64(value);

                        store_[0] = Format.UINT64;
                        *cast(ulong*)&store_[Offset] = temp;
                        stream_.put(store_[0..Offset + ulong.sizeof]);
                    }
                }
            }
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (isSomeChar!T && !is(Unqual!T == enum))
    {
        static if (is(Unqual!T == char)) {
            return pack(cast(ubyte)(value));
        } else static if (is(Unqual!T == wchar)) {
            return pack(cast(ushort)(value));
        } else static if (is(Unqual!T == dchar)) {
            return pack(cast(uint)(value));
        }
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (isFloatingPoint!T && !is(Unqual!T == enum))
    {
        static if (is(Unqual!T == float)) {
            const temp = convertEndianTo!32(_f(value).i);

            store_[0] = Format.FLOAT;
            *cast(uint*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + uint.sizeof]);
        } else static if (is(Unqual!T == double)) {
            const temp = convertEndianTo!64(_d(value).i);

            store_[0] = Format.DOUBLE;
            *cast(ulong*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + ulong.sizeof]);
        } else {
            static if ((real.sizeof > double.sizeof) && EnableReal) {
                store_[0]      = Format.REAL;
                const temp     = _r(value);
                const fraction = convertEndianTo!64(temp.fraction);
                const exponent = convertEndianTo!16(temp.exponent);

                *cast(Unqual!(typeof(fraction))*)&store_[Offset]                   = fraction;
                *cast(Unqual!(typeof(exponent))*)&store_[Offset + fraction.sizeof] = exponent;
                stream_.put(store_[0..$]);
            } else {  // Non-x86 CPUs, real type equals double type.
                pack(cast(double)value);
            }
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (is(Unqual!T == enum))
    {
        pack(cast(OriginalType!T)value);

        return this;
    }


    /// Overload for pack(null) for 2.057 or later
    static if (!is(typeof(null) == void*))
    {
        ref PackerImpl pack(T)(in T value) if (is(Unqual!T == typeof(null)))
        {
            return packNil();
        }
    }


    /// ditto
    ref PackerImpl pack(T)(in T value) if (isPointer!T)
    {
        static if (is(Unqual!T == void*)) {  // for pack(null) for 2.056 or earlier
            enforce(value is null, "Can't serialize void type");
            stream_.put(Format.NIL);
        } else {
            if (value is null)
                stream_.put(Format.NIL);
            else
                pack(mixin(AsteriskOf!T ~ "value"));
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T array) if ((isArray!T || isInstanceOf!(Array, T)) && !is(Unqual!T == enum))
    {
        alias typeof(T.init[0]) U;

        if (array.empty)
            return packNil();

        // Raw bytes
        static if (isByte!(U) || isSomeChar!(U)) {
            ubyte[] raw = cast(ubyte[])array;

            beginRaw(raw.length);
            stream_.put(raw);
        } else {
            beginArray(array.length);
            foreach (elem; array)
                pack(elem);
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(T)(in T array) if (isAssociativeArray!T)
    {
        if (array is null)
            return packNil();

        beginMap(array.length);
        foreach (key, value; array) {
            pack(key);
            pack(value);
        }

        return this;
    }


    /// ditto
    ref PackerImpl pack(Types...)(auto ref const Types objects) if (Types.length > 1)
    {
        foreach (i, T; Types)
            pack(objects[i]);

        return this;
    }


    /**
     * Serializes $(D_PARAM object) and writes to stream.
     *
     * Calling $(D toMsgpack) if $(D_KEYWORD class) and $(D_KEYWORD struct) implement $(D toMsgpack) method. $(D toMsgpack) signature is:
     * -----
     * void toMsgpack(Packer)(ref Packer packer) const
     * -----
     * This method serializes all members of T object if $(D_KEYWORD class) and $(D_KEYWORD struct) don't implement $(D toMsgpack).
     *
     * An object that doesn't implement $(D toMsgpack) is serialized to Array type.
     * -----
     * packer.pack(tuple(true, 1, "Hi!"))  // -> '[true, 1, "Hi!"]', not 'ture, 1, "Hi!"'
     *
     * struct Foo
     * {
     *     int num    = 10;
     *     string msg = "D!";
     * }
     * packer.pack(Foo());  // -> '[10, "D!"]'
     *
     * class Base
     * {
     *     bool flag = true;
     * }
     * class Derived : Base
     * {
     *     double = 0.5f;
     * }
     * packer.pack(new Derived());  // -> '[true, 0.5f]'
     * -----
     *
     * Params:
     *  object = the content to serialize.
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref PackerImpl pack(T)(in T object) if (is(Unqual!T == class))
    {
        if (object is null)
            return packNil();

        static if (hasMember!(T, "toMsgpack"))
        {
            static if (__traits(compiles, { object.toMsgpack(this, withFieldName_); })) {
                object.toMsgpack(this, withFieldName_);
            } else static if (__traits(compiles, { object.toMsgpack(this); })) { // backward compatible
                object.toMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'toMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else {
            if (auto handler = object.classinfo in packHandlers) {
                (*handler)(this, cast(void*)&object);
                return this;
            }
            if (T.classinfo !is object.classinfo) {
                throw new MessagePackException("Can't pack derived class through reference to base class.");
            }

            packObject!(T)(object);
        }

        return this;
    }


    /// ditto
    @trusted
    ref PackerImpl pack(T)(auto ref T object) if (is(Unqual!T == struct) &&
                                                  !isInstanceOf!(Array, T) &&
                                                  !is(Unqual!T == ExtValue))
    {
        static if (hasMember!(T, "toMsgpack"))
        {
            static if (__traits(compiles, { object.toMsgpack(this, withFieldName_); })) {
                object.toMsgpack(this, withFieldName_);
            } else static if (__traits(compiles, { object.toMsgpack(this); })) { // backward compatible
                object.toMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'toMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else static if (isTuple!T) {
            beginArray(object.field.length);
            foreach (f; object.field)
                pack(f);
        } else {  // simple struct
            if (auto handler = typeid(Unqual!T) in packHandlers) {
                (*handler)(this, cast(void*)&object);
                return this;
            }

            immutable memberNum = SerializingMemberNumbers!(T);
            if (withFieldName_)
                beginMap(memberNum);
            else
                beginArray(memberNum);

            if (withFieldName_) {
                foreach (i, f; object.tupleof) {
                    static if (isPackedField!(T.tupleof[i])) {
                        pack(getFieldName!(T, i));
                        static if (hasSerializedAs!(T.tupleof[i])) {
                            alias Proxy = getSerializedAs!(T.tupleof[i]);
                            Proxy.serialize(this, f);
                        } else static if (__traits(compiles, { pack(f); }))
                            pack(f);
                    }
                }
            } else {
                foreach (i, f; object.tupleof) {
                    static if (isPackedField!(T.tupleof[i])) {
                        static if (hasSerializedAs!(T.tupleof[i])) {
                            alias Proxy = getSerializedAs!(T.tupleof[i]);
                            Proxy.serialize(this, f);
                        } else static if (__traits(compiles, { pack(f); }))
                            pack(f);
                    }
                }
            }
        }

        return this;
    }


    void packObject(T)(in T object) if (is(Unqual!T == class))
    {
        alias SerializingClasses!(T) Classes;

        immutable memberNum = SerializingMemberNumbers!(Classes);
        if (withFieldName_)
            beginMap(memberNum);
        else
            beginArray(memberNum);

        foreach (Class; Classes) {
            Class obj = cast(Class)object;
            if (withFieldName_) {
                foreach (i, f ; obj.tupleof) {
                    static if (isPackedField!(Class.tupleof[i])) {
                        pack(getFieldName!(Class, i));
                        static if (hasSerializedAs!(T.tupleof[i])) {
                            alias Proxy = getSerializedAs!(T.tupleof[i]);
                            Proxy.serialize(this, f);
                        } else {
                            pack(f);
                        }
                    }
                }
            } else {
                foreach (i, f ; obj.tupleof) {
                    static if (isPackedField!(Class.tupleof[i])) {
                        static if (hasSerializedAs!(T.tupleof[i])) {
                            alias Proxy = getSerializedAs!(T.tupleof[i]);
                            Proxy.serialize(this, f);
                        } else {
                            pack(f);
                        }
                    }
                }
            }
        }
    }


    /**
     * Serializes the arguments as container to stream.
     *
     * -----
     * packer.packArray(true, 1);  // -> [true, 1]
     * packer.packMap("Hi", 100);  // -> ["Hi":100]
     * -----
     *
     * In packMap, the number of arguments must be even.
     *
     * Params:
     *  objects = the contents to serialize.
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref PackerImpl packArray(Types...)(auto ref const Types objects)
    {
        beginArray(Types.length);
        foreach (i, T; Types)
            pack(objects[i]);
        //pack(objects);  // slow :(

        return this;
    }


    /// ditto
    ref PackerImpl packMap(Types...)(auto ref const Types objects)
    {
        static assert(Types.length % 2 == 0, "The number of arguments must be even");

        beginMap(Types.length / 2);
        foreach (i, T; Types)
            pack(objects[i]);

        return this;
    }

    /**
     * Packs $(D data) as an extended value of $(D type).
     *
     * ----
     * packer.packExt(3, bytes);
     * ----
     *
     * $(D type) must be a signed byte 0-127.
     *
     * Params:
     *  type = the application-defined type for the data
     *  data = an array of bytes
     *
     * Returns:
     *  seld, i.e. for method chaining.
     */
    ref PackerImpl pack(T)(auto ref const T data) if (is(Unqual!T == ExtValue))
    {
        packExt(data.type, data.data);
        return this;
    }

    /**
     * Packs $(D data) as an extended value of $(D type).
     *
     * ----
     * packer.packExt(3, bytes);
     * ----
     *
     * $(D type) must be a signed byte 0-127.
     *
     * Params:
     *  type = the application-defined type for the data
     *  data = an array of bytes
     *
     * Returns:
     *  seld, i.e. for method chaining.
     */
    ref PackerImpl packExt(in byte type, const ubyte[] data) return
    {
        ref PackerImpl packExtFixed(int fmt)
        {
            store_[0] = cast(ubyte)fmt;
            store_[1] = type;
            stream_.put(store_[0 .. 2]);
            stream_.put(data);
            return this;
        }

        // Try packing to a fixed-length type
        if (data.length == 1)
            return packExtFixed(Format.EXT + 0);
        else if (data.length == 2)
            return packExtFixed(Format.EXT + 1);
        else if (data.length == 4)
            return packExtFixed(Format.EXT + 2);
        else if (data.length == 8)
            return packExtFixed(Format.EXT + 3);
        else if (data.length == 16)
            return packExtFixed(Format.EXT + 4);

        int typeByte = void;
        if (data.length <= (2^^8)-1)
        {
            store_[0] = Format.EXT8;
            store_[1] = cast(ubyte)data.length;
            typeByte = 2;

        } else if (data.length <= (2^^16)-1) {
            store_[0] = Format.EXT16;
            const temp = convertEndianTo!16(data.length);
            *cast(ushort*)&store_[Offset] = temp;
            typeByte = 3;
        } else if (data.length <= (2^^32)-1) {
            store_[0] = Format.EXT32;
            const temp = convertEndianTo!32(data.length);
            *cast(uint*)&store_[Offset] = temp;
            typeByte = 5;
        } else
            throw new MessagePackException("Data too large to pack as EXT");

        store_[typeByte] = type;
        stream_.put(store_[0..typeByte+1]);
        stream_.put(data);

        return this;
    }

    /*
     * Serializes raw type-information to stream for binary type.
     */
    void beginRaw(in size_t length)
    {
        import std.conv : text;

        if (length < 32) {
            const ubyte temp = Format.RAW | cast(ubyte)length;
            stream_.put(take8from(temp));
        } else if (length < 65536) {
            const temp = convertEndianTo!16(length);

            store_[0] = Format.RAW16;
            *cast(ushort*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + ushort.sizeof]);
        } else {
            if (length > 0xffffffff)
                throw new MessagePackException(text("size of raw is too long to pack: ", length,  " bytes should be <= ", 0xffffffff));

            const temp = convertEndianTo!32(length);

            store_[0] = Format.RAW32;
            *cast(uint*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + uint.sizeof]);
        }
    }

    /**
     * Serializes the type-information to stream.
     *
     * These methods don't serialize contents.
     * You need to call pack method to serialize contents at your own risk.
     * -----
     * packer.beginArray(3).pack(true, 1);  // -> [true, 1,
     *
     * // other operation
     *
     * packer.pack("Hi!");                  // -> [true, 1, "Hi!"]
     * -----
     *
     * Params:
     *  length = the length of container.
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref PackerImpl beginArray(in size_t length) return
    {
        if (length < 16) {
            const ubyte temp = Format.ARRAY | cast(ubyte)length;
            stream_.put(take8from(temp));
        } else if (length < 65536) {
            const temp = convertEndianTo!16(length);

            store_[0] = Format.ARRAY16;
            *cast(ushort*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + ushort.sizeof]);
        } else {
            const temp = convertEndianTo!32(length);

            store_[0] = Format.ARRAY32;
            *cast(uint*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + uint.sizeof]);
        }

        return this;
    }


    /// ditto
    ref PackerImpl beginMap(in size_t length) return
    {
        if (length < 16) {
            const ubyte temp = Format.MAP | cast(ubyte)length;
            stream_.put(take8from(temp));
        } else if (length < 65536) {
            const temp = convertEndianTo!16(length);

            store_[0] = Format.MAP16;
            *cast(ushort*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + ushort.sizeof]);
        } else {
            const temp = convertEndianTo!32(length);

            store_[0] = Format.MAP32;
            *cast(uint*)&store_[Offset] = temp;
            stream_.put(store_[0..Offset + uint.sizeof]);
        }

        return this;
    }


  private:
    /*
     * Serializes the nil value.
     */
    ref PackerImpl packNil() return
    {
        stream_.put(Format.NIL);
        return this;
    }
}


/// Default serializer
alias PackerImpl!(Appender!(ubyte[])) Packer;  // should be pure struct?


/**
 * Helper for $(D Packer) construction.
 *
 * Params:
 *  stream = the stream to write.
 *  withFieldName = serialize class / struct with field name
 *
 * Returns:
 *  a $(D Packer) object instantiated and initialized according to the arguments.
 */
PackerImpl!(Stream) packer(Stream)(Stream stream, bool withFieldName = false)
{
    return typeof(return)(stream, withFieldName);
}


version(unittest)
{
    package import std.file, core.stdc.string;

    package mixin template DefinePacker()
    {
        Packer packer;
    }

    package mixin template DefineDictionalPacker()
    {
        Packer packer = Packer(false);
    }
}


unittest
{
    { // unique value
        mixin DefinePacker;

        ubyte[] result = [Format.NIL, Format.TRUE, Format.FALSE];

        packer.pack(null, true, false);
        foreach (i, value; packer.stream.data)
            assert(value == result[i]);
    }
    { // uint *
        static struct UTest { ubyte format; ulong value; }

        enum : ulong { A = ubyte.max, B = ushort.max, C = uint.max, D = ulong.max }

        static UTest[][] utests = [
            [{Format.UINT8, A}],
            [{Format.UINT8, A}, {Format.UINT16, B}],
            [{Format.UINT8, A}, {Format.UINT16, B}, {Format.UINT32, C}],
            [{Format.UINT8, A}, {Format.UINT16, B}, {Format.UINT32, C}, {Format.UINT64, D}],
        ];

        foreach (I, T; TypeTuple!(ubyte, ushort, uint, ulong)) {
            foreach (i, test; utests[I]) {
                mixin DefinePacker;

                packer.pack(cast(T)test.value);
                assert(packer.stream.data[0] == test.format);

                switch (i) {
                case 0:
                    auto answer = take8from!(T.sizeof * 8)(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, ubyte.sizeof) == 0);
                    break;
                case 1:
                    auto answer = convertEndianTo!16(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, ushort.sizeof) == 0);
                    break;
                case 2:
                    auto answer = convertEndianTo!32(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, uint.sizeof) == 0);
                    break;
                default:
                    auto answer = convertEndianTo!64(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, ulong.sizeof) == 0);
                }
            }
        }
    }
    { // int *
        static struct STest { ubyte format; long value; }

        enum : long { A = byte.min, B = short.min, C = int.min, D = long.min }

        static STest[][] stests = [
            [{Format.INT8, A}],
            [{Format.INT8, A}, {Format.INT16, B}],
            [{Format.INT8, A}, {Format.INT16, B}, {Format.INT32, C}],
            [{Format.INT8, A}, {Format.INT16, B}, {Format.INT32, C}, {Format.INT64, D}],
        ];

        foreach (I, T; TypeTuple!(byte, short, int, long)) {
            foreach (i, test; stests[I]) {
                mixin DefinePacker;

                packer.pack(cast(T)test.value);
                assert(packer.stream.data[0] == test.format);

                switch (i) {
                case 0:
                    auto answer = take8from!(T.sizeof * 8)(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, byte.sizeof) == 0);
                    break;
                case 1:
                    auto answer = convertEndianTo!16(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, short.sizeof) == 0);
                    break;
                case 2:
                    auto answer = convertEndianTo!32(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, int.sizeof) == 0);
                    break;
                default:
                    auto answer = convertEndianTo!64(test.value);
                    assert(memcmp(&packer.stream.data[1], &answer, long.sizeof) == 0);
                }
            }
        }
    }
    { // fload, double
        static if ((real.sizeof == double.sizeof) || !EnableReal)
        {
            alias TypeTuple!(float, double, double) FloatingTypes;
            static struct FTest { ubyte format; double value; }

            static FTest[] ftests = [
                {Format.FLOAT,  float.min_normal},
                {Format.DOUBLE, double.max},
                {Format.DOUBLE, double.max},
            ];
        }
        else
        {
            alias TypeTuple!(float, double, real) FloatingTypes;
            static struct FTest { ubyte format; real value; }

            static FTest[] ftests = [
                {Format.FLOAT,  float.min_normal},
                {Format.DOUBLE, double.max},
                {Format.REAL,   real.max},
            ];
        }

        foreach (I, T; FloatingTypes) {
            mixin DefinePacker;

            packer.pack(cast(T)ftests[I].value);
            assert(packer.stream.data[0] == ftests[I].format);

            switch (I) {
            case 0:
                const answer = convertEndianTo!32(_f(cast(T)ftests[I].value).i);
                assert(memcmp(&packer.stream.data[1], &answer, float.sizeof) == 0);
                break;
            case 1:
                const answer = convertEndianTo!64(_d(cast(T)ftests[I].value).i);
                assert(memcmp(&packer.stream.data[1], &answer, double.sizeof) == 0);
                break;
            default:
                static if (EnableReal)
                {
                    const t = _r(cast(T)ftests[I].value);
                    const f = convertEndianTo!64(t.fraction);
                    const e = convertEndianTo!16(t.exponent);
                    assert(memcmp(&packer.stream.data[1],            &f, f.sizeof) == 0);
                    assert(memcmp(&packer.stream.data[1 + f.sizeof], &e, e.sizeof) == 0);
                }
                else
                {
                    const answer = convertEndianTo!64(_d(cast(T)ftests[I].value).i);
                    assert(memcmp(&packer.stream.data[1], &answer, double.sizeof) == 0);
                }
            }
        }
    }
    { // pointer
        static struct PTest
        {
            ubyte format;

            union
            {
                ulong*  p0;
                long*   p1;
                double* p2;
            }
        }

        PTest[] ptests = [PTest(Format.UINT64), PTest(Format.INT64), PTest(Format.DOUBLE)];

        ulong  v0 = ulong.max;
        long   v1 = long.min;
        double v2 = double.max;

        foreach (I, Index; TypeTuple!("0", "1", "2")) {
            mixin DefinePacker;

            mixin("ptests[I].p" ~ Index ~ " = &v" ~ Index ~ ";");

            packer.pack(mixin("ptests[I].p" ~ Index));
            assert(packer.stream.data[0] == ptests[I].format);

            switch (I) {
            case 0:
                auto answer = convertEndianTo!64(*ptests[I].p0);
                assert(memcmp(&packer.stream.data[1], &answer, ulong.sizeof) == 0);
                break;
            case 1:
                auto answer = convertEndianTo!64(*ptests[I].p1);
                assert(memcmp(&packer.stream.data[1], &answer, long.sizeof) == 0);
                break;
            default:
                const answer = convertEndianTo!64(_d(*ptests[I].p2).i);
                assert(memcmp(&packer.stream.data[1], &answer, double.sizeof) == 0);
            }
        }
    }
    { // enum
        enum E : ubyte { A = ubyte.max }

        mixin DefinePacker; E e = E.A;

        packer.pack(e);
        assert(packer.stream.data[0] == Format.UINT8);

        auto answer = E.A;
        assert(memcmp(&packer.stream.data[1], &answer, (OriginalType!E).sizeof) == 0);
    }
    { // enum with string
        enum E2 : string { A = "test" }

        mixin DefinePacker; E2 e = E2.A;

        packer.pack(e);
        assert(packer.stream.data[0] == (Format.RAW | 0x04));
    }
    { // container
        static struct CTest { ubyte format; size_t value; }

        enum : ulong { A = 16 / 2, B = ushort.max, C = uint.max }

        static CTest[][] ctests = [
            [{Format.ARRAY | A, Format.ARRAY | A}, {Format.ARRAY16, B}, {Format.ARRAY32, C}],
            [{Format.MAP   | A, Format.MAP   | A}, {Format.MAP16,   B}, {Format.MAP32,   C}],
            [{Format.RAW   | A, Format.RAW   | A}, {Format.RAW16,   B}, {Format.RAW32,   C}],
        ];

        foreach (I, Name; TypeTuple!("Array", "Map", "Raw")) {
            auto test = ctests[I];

            foreach (i, T; TypeTuple!(ubyte, ushort, uint)) {
                mixin DefinePacker;
                mixin("packer.begin" ~ Name ~ "(i ? test[i].value : A);");

                assert(packer.stream.data[0] == test[i].format);

                switch (i) {
                case 0:
                    auto answer = take8from(test[i].value);
                    assert(memcmp(&packer.stream.data[0], &answer, ubyte.sizeof) == 0);
                    break;
                case 1:
                    auto answer = convertEndianTo!16(test[i].value);
                    assert(memcmp(&packer.stream.data[1], &answer, ushort.sizeof) == 0);
                    break;
                default:
                    auto answer = convertEndianTo!32(test[i].value);
                    assert(memcmp(&packer.stream.data[1], &answer, uint.sizeof) == 0);
                }
            }
        }
    }

    version (X86_64) // can't create a long enough array to trigger this on x86
    { // larger spec size for string / binary
        mixin DefinePacker;

        try {
            // using malloc because - hopefully - this means we don't
            // actually physically allocate such a huge amount of memory
            import core.stdc.stdlib;
            auto len = 0xffffffffUL + 1;
            auto bins = (cast(byte*)malloc(len))[0 .. len];
            assert(bins);
            scope(exit) free(bins.ptr);
            packer.pack(bins);
            assert(false); //check it wasn't allowed
        } catch (MessagePackException e) {
        }
    }
    { // user defined
        {
            static struct S
            {
                uint num = uint.max;

                void toMsgpack(P)(ref P p) const { p.packArray(num); }
            }

            mixin DefinePacker; S test;

            packer.pack(test);

            assert(packer.stream.data[0] == (Format.ARRAY | 1));
            assert(packer.stream.data[1] ==  Format.UINT32);
            assert(memcmp(&packer.stream.data[2], &test.num, uint.sizeof) == 0);
        }
        {
            mixin DefinePacker; auto test = tuple(true, false, uint.max);

            packer.pack(test);

            assert(packer.stream.data[0] == (Format.ARRAY | 3));
            assert(packer.stream.data[1] ==  Format.TRUE);
            assert(packer.stream.data[2] ==  Format.FALSE);
            assert(packer.stream.data[3] ==  Format.UINT32);
            assert(memcmp(&packer.stream.data[4], &test.field[2], uint.sizeof) == 0);
        }
        {
            static class C
            {
                uint num;

                this(uint n) { num = n; }

                void toMsgpack(P)(ref P p) const { p.packArray(num); }
            }

            mixin DefinePacker; C test = new C(ushort.max);

            packer.pack(test);

            assert(packer.stream.data[0] == (Format.ARRAY | 1));
            assert(packer.stream.data[1] ==  Format.UINT16);
            assert(memcmp(&packer.stream.data[2], &test.num, ushort.sizeof) == 0);
        }
    }
    { // simple struct and class
        {
            static struct Simple
            {
                uint num = uint.max;
            }

            static struct SimpleWithNonPacked1
            {
                uint num = uint.max;
                @nonPacked string str = "ignored";
            }

            static struct SimpleWithNonPacked2
            {
                @nonPacked string str = "ignored";
                uint num = uint.max;
            }

            static struct SimpleWithSkippedTypes
            {
                int function(int) fn;
                int delegate(int) dg;
                uint num = uint.max;
            }

            foreach (Type; TypeTuple!(Simple, SimpleWithNonPacked1, SimpleWithNonPacked2, SimpleWithSkippedTypes)) {
                mixin DefinePacker;

                Type test;
                packer.pack(test);

                assert(packer.stream.data[0] == (Format.ARRAY | 1));
                assert(packer.stream.data[1] ==  Format.UINT32);
                assert(memcmp(&packer.stream.data[2], &test.num, uint.sizeof) == 0);
            }
        }

        static class SimpleA
        {
            bool flag = true;
        }

        static class SimpleB : SimpleA
        {
            ubyte type = 100;
        }

        static class SimpleC : SimpleB
        {
            uint num = uint.max;
        }

        static class SimpleCWithNonPacked1 : SimpleB
        {
            uint num = uint.max;
            @nonPacked string str = "ignored";
        }

        static class SimpleCWithNonPacked2 : SimpleB
        {
            @nonPacked string str = "ignored";
            uint num = uint.max;
        }

        static class SimpleCWithSkippedTypes : SimpleB
        {
            uint num = uint.max;
            int function(int) fn;
            int delegate(int) dg;
        }

        {  // from derived class
            foreach (Type; TypeTuple!(SimpleC, SimpleCWithNonPacked1, SimpleCWithNonPacked2, SimpleCWithSkippedTypes)) {
                mixin DefinePacker;

                Type test = new Type();
                packer.pack(test);

                assert(packer.stream.data[0] == (Format.ARRAY | 3));
                assert(packer.stream.data[1] ==  Format.TRUE);
                assert(packer.stream.data[2] ==  100);
                assert(packer.stream.data[3] ==  Format.UINT32);
                assert(memcmp(&packer.stream.data[4], &test.num, uint.sizeof) == 0);
            }
        }
        {  // from base class
            mixin DefinePacker; SimpleB test = new SimpleC();

            try {
                packer.pack(test);
                assert(false);
            } catch (Exception e) { }
        }
    }

    // ext types
    {
        byte type = 7; // an arbitrary type value

        // fixexts
        {
            ubyte[16] data;
            data[] = 1;
            foreach (L; TypeTuple!(1, 2, 4, 8, 16))
            {
                mixin DefinePacker;
                packer.pack(ExtValue(type, data[0 .. L]));

                // format, type, data
                assert(packer.stream.data.length == 2 + L);
                const l = 2 ^^ (packer.stream.data[0] - Format.EXT);
                assert(l == L);
                assert(packer.stream.data[1] == type);
                assert(packer.stream.data[2 .. 2+l] == data[0 .. L]);
            }
        }

        // ext8
        {
            foreach (L; TypeTuple!(3, 7, 255))
            {
                ubyte[] data = new ubyte[](L);
                data[] = 1;

                mixin DefinePacker;
                packer.pack(ExtValue(type, data[0 .. L]));

                // format, length, type, data
                assert(packer.stream.data.length == 3 + L);
                assert(packer.stream.data[0] == Format.EXT8);
                assert(packer.stream.data[1] == L);
                assert(packer.stream.data[2] == type);
                assert(packer.stream.data[3 .. 3 + L] == data);
            }
        }

        // ext16
        {
            foreach (L; TypeTuple!(256, (2^^16)-1))
            {
                ubyte[] data = new ubyte[](L);
                data[] = 1;

                mixin DefinePacker;
                packer.pack(ExtValue(type, data[0 .. L]));

                // format, length, type, data
                import std.conv : text;
                assert(packer.stream.data.length == 4 + L, text(packer.stream.data.length));
                assert(packer.stream.data[0] == Format.EXT16);

                ushort l = convertEndianTo!16(L);
                assert(memcmp(&packer.stream.data[1], &l, ushort.sizeof) == 0);
                assert(packer.stream.data[3] == type);
                assert(packer.stream.data[4 .. 4 + L] == data);
            }
        }

        // ext32
        {
            foreach (L; TypeTuple!(2^^16, 2^^17))
            {
                ubyte[] data = new ubyte[](L);
                data[] = 1;

                mixin DefinePacker;
                packer.pack(ExtValue(type, data[0 .. L]));

                // format, length, type, data
                import std.conv : text;
                assert(packer.stream.data.length == 6 + L, text(packer.stream.data.length));
                assert(packer.stream.data[0] == Format.EXT32);

                uint l = convertEndianTo!32(L);
                assert(memcmp(&packer.stream.data[1], &l, uint.sizeof) == 0);
                assert(packer.stream.data[5] == type);
                assert(packer.stream.data[6 .. 6 + L] == data);
            }
        }
    }
}
