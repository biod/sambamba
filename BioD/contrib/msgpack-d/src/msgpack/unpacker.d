module msgpack.unpacker;

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


// for unpack without calling constructor
private extern(C) Object _d_newclass(const ClassInfo);


/**
 * This $(D Unpacker) is a $(D MessagePack) direct-conversion deserializer
 *
 * This implementation is suitable for fixed data.
 *
 * Example:
 * -----
 * // serializedData is [10, 0.1, false]
 * auto unpacker = Unpacker(serializedData);
 *
 * uint   n;
 * double d;
 * bool   b;
 *
 * unpacker.unpackArray(n, d, b);
 *
 * // using Tuple
 * Tuple!(uint, double, bool) record;
 * unpacker.unpack(record);  // record is [10, 0.1, false]
 * -----
 *
 * NOTE:
 *  Unpacker becomes template struct if Phobos supports truly IO module.
 */
struct Unpacker
{
  private:
    static @system
    {
        alias void delegate(ref Unpacker, void*) UnpackHandler;
        UnpackHandler[TypeInfo] unpackHandlers;

        public void registerHandler(T, alias Handler)()
        {
            unpackHandlers[typeid(T)] = delegate(ref Unpacker unpacker, void* obj) {
                Handler(unpacker, *cast(T*)obj);
            };
        }

        public void register(T)()
        {
            unpackHandlers[typeid(T)] = delegate(ref Unpacker unpacker, void* obj) {
                unpacker.unpackObject(*cast(T*)obj);
            };
        }

    }

    enum Offset = 1;

    mixin InternalBuffer;

    bool withFieldName_;


  public:
    /**
     * Constructs a $(D Unpacker).
     *
     * Params:
     *  target     = byte buffer to deserialize
     *  bufferSize = size limit of buffer size
     */
    this(in ubyte[] target, in size_t bufferSize = 8192, bool withFieldName = false)
    {
        initializeBuffer(target, bufferSize);
        withFieldName_ = withFieldName;
    }


    /**
     * Clears states for next deserialization.
     */
    @safe
    nothrow void clear()
    {
        used_ = offset_ = parsed_ = 0;
        hasRaw_ = false;
    }


    /**
     * Deserializes $(D_PARAM T) object and assigns to $(D_PARAM value).
     *
     * If the argument is pointer, dereferences pointer and assigns deserialized value.
     * -----
     * int* a;
     * unpacker.unpack(a)  // enforce throws Exception because a is null or
     *                     // no throw if deserialized value is nil
     *
     * int b; a = &b;
     * unpacker.unpack(b)  // b is deserialized value or
     *                     // assigns null if deserialized value is nil
     * -----
     *
     * Params:
     *  value = the reference of value to assign.
     *
     * Returns:
     *  self, i.e. for method chaining.
     *
     * Throws:
     *  UnpackException when doesn't read from buffer or precision loss occurs and
     *  MessagePackException when $(D_PARAM T) type doesn't match serialized type.
     */
    ref Unpacker unpack(T)(ref T value) if (is(Unqual!T == bool))
    {
        canRead(Offset, 0);
        const header = read();

        switch (header) {
        case Format.TRUE:
            value = true;
            break;
        case Format.FALSE:
            value = false;
            break;
        default:
            rollback(0, "bool", cast(Format)header);
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (isUnsigned!T && !is(Unqual!T == enum))
    {
        canRead(Offset, 0);
        const header = read();

        if (0x00 <= header && header <= 0x7f) {
            value = header;
        } else {
            switch (header) {
            case Format.UINT8:
                canRead(ubyte.sizeof);
                value = read();
                break;
            case Format.UINT16:
                canRead(ushort.sizeof);
                auto us = load16To!ushort(read(ushort.sizeof));
                if (us > T.max)
                    rollback(ushort.sizeof, T.stringof, Format.UINT16);
                value = cast(T)us;
                break;
            case Format.UINT32:
                canRead(uint.sizeof);
                auto ui = load32To!uint(read(uint.sizeof));
                if (ui > T.max)
                    rollback(uint.sizeof, T.stringof, Format.UINT32);
                value = cast(T)ui;
                break;
            case Format.UINT64:
                canRead(ulong.sizeof);
                auto ul = load64To!ulong(read(ulong.sizeof));
                if (ul > T.max)
                    rollback(ulong.sizeof, T.stringof, Format.UINT64);
                value = cast(T)ul;
                break;
            default:
                rollback(0, T.stringof, cast(Format)header);
            }
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (isSigned!T && isIntegral!T && !is(Unqual!T == enum))
    {
        canRead(Offset, 0);
        const header = read();

        if (0x00 <= header && header <= 0x7f) {
            value = cast(T)header;
        } else if (0xe0 <= header && header <= 0xff) {
            auto b = cast(T)(cast(ubyte)(-int(header)));
            static if(T.sizeof < int.sizeof)
                value = cast(T)(-int(b));
            else
                value = -b;
        } else {
            switch (header) {
            case Format.UINT8:
                canRead(ubyte.sizeof);
                auto ub = read();
                if (ub > T.max)
                    rollback(ubyte.sizeof, T.stringof, Format.UINT8);
                value = cast(T)ub;
                break;
            case Format.UINT16:
                canRead(ushort.sizeof);
                auto us = load16To!ushort(read(ushort.sizeof));
                if (us > T.max)
                    rollback(ushort.sizeof, T.stringof, Format.UINT16);
                value = cast(T)us;
                break;
            case Format.UINT32:
                canRead(uint.sizeof);
                auto ui = load32To!uint(read(uint.sizeof));
                if (ui > T.max)
                    rollback(uint.sizeof, T.stringof, Format.UINT32);
                value = cast(T)ui;
                break;
            case Format.UINT64:
                canRead(ulong.sizeof);
                auto ul = load64To!ulong(read(ulong.sizeof));
                if (ul > T.max)
                    rollback(ulong.sizeof, T.stringof, Format.UINT64);
                value = cast(T)ul;
                break;
            case Format.INT8:
                canRead(byte.sizeof);
                value = cast(byte)read();
                break;
            case Format.INT16:
                canRead(short.sizeof);
                auto s = load16To!short(read(short.sizeof));
                if (s < T.min || T.max < s)
                    rollback(short.sizeof, T.stringof, Format.INT16);
                value = cast(T)s;
                break;
            case Format.INT32:
                canRead(int.sizeof);
                auto i = load32To!int(read(int.sizeof));
                if (i < T.min || T.max < i)
                    rollback(int.sizeof, T.stringof, Format.INT32);
                value = cast(T)i;
                break;
            case Format.INT64:
                canRead(long.sizeof);
                auto l = load64To!long(read(long.sizeof));
                if (l < T.min || T.max < l)
                    rollback(long.sizeof, T.stringof, Format.INT64);
                value = cast(T)l;
                break;
            default:
                rollback(0, T.stringof, cast(Format)header);
            }
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (isSomeChar!T && !is(Unqual!T == enum))
    {
        static if (is(Unqual!T == char)) {
            ubyte tmp;
        } else static if (is(Unqual!T == wchar)) {
            ushort tmp;
        } else static if (is(Unqual!T == dchar)) {
            uint tmp;
        }
        unpack(tmp);
        value = cast(T)(tmp);
        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (isFloatingPoint!T && !is(Unqual!T == enum))
    {
        canRead(Offset, 0);
        const header = read();

        switch (header) {
        case Format.FLOAT:
            _f temp;

            canRead(uint.sizeof);
            temp.i = load32To!uint(read(uint.sizeof));
            value  = temp.f;
            break;
        case Format.DOUBLE:
            // check precision loss
            static if (is(Unqual!T == float))
                rollback(0, T.stringof, Format.DOUBLE);

            _d temp;

            canRead(ulong.sizeof);
            temp.i = load64To!ulong(read(ulong.sizeof));
            value  = temp.f;
            break;
        case Format.REAL:
            static if (!EnableReal)
            {
                rollback(0, "real is disabled", Format.REAL);
            }
            else
            {
                // check precision loss
                static if (is(Unqual!T == float) || is(Unqual!T == double))
                    rollback(0, T.stringof, Format.REAL);

                canRead(RealSize);

                version (NonX86)
                {
                    CustomFloat!80 temp;

                    const frac = load64To!ulong (read(ulong.sizeof));
                    const exp  = load16To!ushort(read(ushort.sizeof));

                    temp.significand = frac;
                    temp.exponent    = exp & 0x7fff;
                    temp.sign        = exp & 0x8000 ? true : false;

                    // NOTE: temp.get!real is inf on non-x86 when deserialized value is larger than double.max.
                    value = temp.get!real;
                }
                else
                {
                    _r temp;

                    temp.fraction = load64To!(typeof(temp.fraction))(read(temp.fraction.sizeof));
                    temp.exponent = load16To!(typeof(temp.exponent))(read(temp.exponent.sizeof));

                    value = temp.f;
                }
            }

            break;
        default:
            rollback(0, T.stringof, cast(Format)header);
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (is(Unqual!T == enum))
    {
        OriginalType!T temp;

        unpack(temp);

        value = cast(T)temp;

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(T value) if (isPointer!T)
    {
        static if (is(Unqual!T == void*)) {
            enforce(value !is null,  "Can't deserialize void type");
            unpackNil(value);
        } else {
            if (checkNil())
                unpackNil(value);
            else
                enforce(value !is null, T.stringof ~ " is null pointer");

            unpack(mixin(AsteriskOf!T ~ "value"));
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T value) if (is(T == ExtValue))
    {
        canRead(Offset, 0);
        const header = read();

        // Fixed
        if (header >= Format.EXT && header <= Format.EXT + 4)
        {
            const length = 2^^(header - Format.EXT);
            canRead(1 + length);

            value.type = read();
            value.data = read(length);
            return this;
        }

        // Dynamic length
        uint length;
        switch (header) with (Format)
        {
            case EXT8:
                canRead(1);
                length = read();
                break;
            case EXT16:
                canRead(2);
                length = load16To!ushort(read(2));
                break;
            case EXT32:
                canRead(4);
                length = load32To!uint(read(4));
                break;
            default:
                rollback(0, T.stringof, cast(Format)header);
        }

        canRead(1 + length);
        value.type = read();
        value.data = read(length);

        return this;
    }


    /// ditto
    ref Unpacker unpack(Types...)(ref Types objects) if (Types.length > 1)
    {
        foreach (i, T; Types)
            unpack!(T)(objects[i]);

        return this;
    }


    /**
     * Deserializes $(D_PARAM T) object and assigns to $(D_PARAM array).
     *
     * This is convenient method for array deserialization.
     * Rollback will be completely successful if you deserialize raw type((u)byte[] or string types).
     * But, Rollback will be one element(e.g. int) if you deserialize other types(e.g. int[], int[int])
     *
     * No assign if the length of deserialized object is 0.
     *
     * In a static array, this method checks the length. Do rollback and throw exception
     * if length of $(D_PARAM array) is different from length of deserialized object.
     *
     * Params:
     *  array = the reference of array to assign.
     *
     * Returns:
     *  self, i.e. for method chaining.
     *
     * Throws:
     *  UnpackException when doesn't read from buffer or precision loss occurs and
     *  MessagePackException when $(D_PARAM T) type doesn't match serialized type.
     */
    ref Unpacker unpack(T)(ref T array) if ((isArray!T ||
                                             isInstanceOf!(Array, T)) &&
                                            !is(Unqual!T == enum))
    {
        alias typeof(T.init[0]) U;

        /*
         * Deserializes type-information of raw type.
         */
        @safe
        size_t beginRaw()
        {
            canRead(Offset, 0);
            const  header = read();
            size_t length;

            if (0xa0 <= header && header <= 0xbf) {
                length = header & 0x1f;
            } else {
                switch (header) {
                case Format.BIN8, Format.STR8:
                    canRead(ubyte.sizeof);
                    length = read();
                    break;
                case Format.BIN16, Format.RAW16:
                    canRead(ushort.sizeof);
                    length = load16To!size_t(read(ushort.sizeof));
                    break;
                case Format.BIN32, Format.RAW32:
                    canRead(uint.sizeof);
                    length = load32To!size_t(read(uint.sizeof));
                    break;
                case Format.NIL:
                    break;
                default:
                    rollback(0, T.stringof, cast(Format)header);
                }
            }

            return length;
        }


        if (checkNil()) {
            static if (isStaticArray!T) {
                onInvalidType("static array", Format.NIL);
            } else {
                return unpackNil(array);
            }
        }

        // Raw bytes
        static if (isByte!U || isSomeChar!U)
            auto length = beginRaw();
        else
            auto length = beginArray();

        if(length > buffer_.length) {
            import std.conv: text;
            throw new MessagePackException(text("Invalid array size in byte stream: Length (", length,
                                                ") is larger than internal buffer size (", buffer_.length, ")"));
        }

        // Raw bytes
        static if (isByte!U || isSomeChar!U) {
            auto offset = calculateSize!(true)(length);
            if (length == 0)
                return this;

            static if (isStaticArray!T) {
                if (length != array.length)
                    rollback(offset, "static array was given but the length is mismatched");
            }

            canRead(length, offset + Offset);
            static if (isStaticArray!T) {
                array[] = (cast(U[])read(length))[0 .. T.length];
            } else {
                array = cast(T)read(length);
            }

            static if (isDynamicArray!T)
                hasRaw_ = true;
        } else {
            if (length == 0)
                return this;

            static if (isStaticArray!T) {
                if (length != array.length)
                    rollback(calculateSize(length), "static array was given but the length is mismatched");
            } else {
                array.length = length;
            }

            foreach (i; 0..length)
                unpack(array[i]);
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T array) if (isAssociativeArray!T)
    {
        alias typeof(T.init.keys[0])   K;
        alias typeof(T.init.values[0]) V;

        if (checkNil())
            return unpackNil(array);

        auto length = beginMap();
        if (length == 0)
            return this;

        foreach (i; 0..length) {
            K k; unpack(k);
            V v; unpack(v);
            array[k] = v;
        }

        return this;
    }

    /**
     * Deserializes $(D_PARAM T) object and assigns to $(D_PARAM object).
     *
     * Calling $(D fromMsgpack) if $(D_KEYWORD class) and $(D_KEYWORD struct) implement $(D fromMsgpack) method. $(D fromMsgpack) signature is:
     * -----
     * void fromMsgpack(ref Unpacker unpacker)
     * -----
     * Assumes $(D std.typecons.Tuple) or simple struct if $(D_KEYWORD struct) doesn't implement $(D fromMsgpack).
     * Checks length if $(D_PARAM T) is a $(D std.typecons.Tuple) or simple struct.
     *
     * Params:
     *  object = the reference of object to assign.
     *  args   = the arguments to class constructor(class only).
     *           This is used at new statement if $(D_PARAM object) is $(D_KEYWORD null).
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref Unpacker unpack(T, Args...)(ref T object, auto ref Args args) if (is(Unqual!T == class))
    {
        if (checkNil())
            return unpackNil(object);

        if (object is null) {
            static if (Args.length == 0) {
                static if (__traits(compiles, { new T(); }))
                    object = new T();
                else
                    object = cast(T)_d_newclass(T.classinfo);
            } else static if (__traits(compiles, { new T(args); })) {
                object = new T(args);
            } else {
                throw new MessagePackException("Don't know how to construct class type '" ~ Unqual!T.stringof ~ "' with argument types '" ~ Args.stringof ~ "'.");
            }
        }

        static if (hasMember!(T, "fromMsgpack"))
        {
            static if (__traits(compiles, { object.fromMsgpack(this, withFieldName_); })) {
              object.fromMsgpack(this, withFieldName_);
            } else static if (__traits(compiles, { object.fromMsgpack(this); })) { // backward compatible
                object.fromMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'fromMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else {
            if (auto handler = object.classinfo in unpackHandlers) {
                (*handler)(this, cast(void*)&object);
                return this;
            }
            if (T.classinfo !is object.classinfo) {
                throw new MessagePackException("Can't unpack derived class through reference to base class.");
            }

            unpackObject(object);
        }

        return this;
    }


    /// ditto
    ref Unpacker unpack(T)(ref T object) if (is(Unqual!T == struct) &&
                                             !is(Unqual!T == ExtValue))
    {
        static if (hasMember!(T, "fromMsgpack"))
        {
            static if (__traits(compiles, { object.fromMsgpack(this); })) {
                object.fromMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'fromMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else {
            if (auto handler = typeid(T) in unpackHandlers) {
                (*handler)(this, cast(void*)&object);
                return this;
            }

            size_t length = withFieldName_ ? beginMap() : beginArray();
            if (length == 0)
                return this;

            static if (isTuple!T) {
                if (length != T.Types.length)
                    rollback(calculateSize(length), "the number of tuple fields is mismatched");

                foreach (i, Type; T.Types)
                    unpack(object.field[i]);
            } else {  // simple struct
                //if (length != object.tupleof.length)
                if (length != SerializingMemberNumbers!(T))
                    rollback(calculateSize(length), "the number of struct fields is mismatched");

                if (withFieldName_) {
                    foreach (i, member; object.tupleof) {
                        static if (isPackedField!(T.tupleof[i]))
                        {
                            string fieldName;
                            unpack(fieldName);

                            if (fieldName == getFieldName!(T, i)) {
                                static if (hasSerializedAs!(T.tupleof[i])) {
                                    alias Proxy = getSerializedAs!(T.tupleof[i]);
                                    Proxy.deserialize(this, object.tupleof[i]);
                                } else {
                                    unpack(object.tupleof[i]);
                                }
                            } else {
                                assert(false, "Invalid field name: '" ~ fieldName ~ "', expect '" ~ getFieldName!(T, i) ~ "'");
                            }
                        }
                    }
                } else {
                    foreach (i, member; object.tupleof) {
                        static if (isPackedField!(T.tupleof[i])) {
                            static if (hasSerializedAs!(T.tupleof[i])) {
                                alias Proxy = getSerializedAs!(T.tupleof[i]);
                                Proxy.deserialize(this, object.tupleof[i]);
                            } else {
                                unpack(object.tupleof[i]);
                            }
                        }
                    }
                }
            }
        }

        return this;
    }


    void unpackObject(T)(ref T object) if (is(Unqual!T == class))
    {
        alias SerializingClasses!(T) Classes;

        size_t length = withFieldName_ ? beginMap() : beginArray();
        if (length == 0)
            return;

        if (length != SerializingMemberNumbers!(Classes))
            rollback(calculateSize(length),  "the number of class fields is mismatched");

        if (withFieldName_) {
            foreach (_; 0..length) {
                string fieldName;
                unpack(fieldName);

                foreach (Class; Classes) {
                    Class obj = cast(Class)object;

                    foreach (i, member; obj.tupleof) {
                        static if (isPackedField!(Class.tupleof[i]))
                        {
                            if (fieldName == getFieldName!(Class, i)) {
                                static if (hasSerializedAs!(Class.tupleof[i])) {
                                    alias Proxy = getSerializedAs!(Class.tupleof[i]);
                                    Proxy.deserialize(this, obj.tupleof[i]);
                                } else {
                                    unpack(obj.tupleof[i]);
                                }
                                goto endLoop;
                            }
                        }
                    }
                }
                assert(false, "Invalid field name: '" ~ fieldName~"' ");

            endLoop:
                continue;
            }
        } else {
            foreach (Class; Classes) {
                Class obj = cast(Class)object;

                foreach (i, member; obj.tupleof) {
                    static if (isPackedField!(Class.tupleof[i])) {
                        static if (hasSerializedAs!(Class.tupleof[i])) {
                            alias Proxy = getSerializedAs!(Class.tupleof[i]);
                            Proxy.deserialize(this, obj.tupleof[i]);
                        } else {
                            unpack(obj.tupleof[i]);
                        }
                    }
                }
            }
        }
    }


    /**
     * Deserializes the container object and assigns to each argument.
     *
     * These methods check the length. Do rollback if
     * the length of arguments is different from length of deserialized object.
     *
     * In unpackMap, the number of arguments must be even.
     *
     * Params:
     *  objects = the references of object to assign.
     *
     * Returns:
     *  self, i.e. for method chaining.
     */
    ref Unpacker unpackArray(Types...)(ref Types objects)
    {
        auto length = beginArray();
        if (length != Types.length)
            rollback(calculateSize(length), "the number of deserialized objects is mismatched");

        foreach (i, T; Types)
            unpack(objects[i]);
        // unpack(objects);  // slow :(

        return this;
    }


    /// ditto
    ref Unpacker unpackMap(Types...)(ref Types objects)
    {
        static assert(Types.length % 2 == 0, "The number of arguments must be even");

        auto length = beginMap();
        if (length != Types.length / 2)
            rollback(calculateSize(length), "the number of deserialized objects is mismatched");

        foreach (i, T; Types)
            unpack(objects[i]);

        return this;
    }


    /**
     * Deserializes the type-information of container.
     *
     * These methods don't deserialize contents.
     * You need to call unpack method to deserialize contents at your own risk.
     * -----
     * // serialized data is [1, "Hi!"];
     * int num;
     * unpacker.beginArray(2).unpack(num);  // num is 1
     *
     * // other operation
     *
     * string str;
     * unpacker.unpack(str);  // str is "Hi!"
     * -----
     *
     * Returns:
     *  the container size.
     */
    @safe
    size_t beginArray()
    {
        canRead(Offset, 0);
        const  header = read();
        size_t length;

        if (0x90 <= header && header <= 0x9f) {
            length = header & 0x0f;
        } else {
            switch (header) {
            case Format.ARRAY16:
                canRead(ushort.sizeof);
                length = load16To!size_t(read(ushort.sizeof));
                break;
            case Format.ARRAY32:
                canRead(uint.sizeof);
                length = load32To!size_t(read(uint.sizeof));
                break;
            case Format.NIL:
                break;
            default:
                rollback(0, "array", cast(Format)header);
            }
        }

        return length;
    }


    /// ditto
    @safe
    size_t beginMap()
    {
        canRead(Offset, 0);
        const  header = read();
        size_t length;

        if (0x80 <= header && header <= 0x8f) {
            length = header & 0x0f;
        } else {
            switch (header) {
            case Format.MAP16:
                canRead(ushort.sizeof);
                length = load16To!size_t(read(ushort.sizeof));
                break;
            case Format.MAP32:
                canRead(uint.sizeof);
                length = load32To!size_t(read(uint.sizeof));
                break;
            case Format.NIL:
                break;
            default:
                rollback(0, "map", cast(Format)header);
            }
        }

        return length;
    }


    /**
     * Unpacks an EXT value into $(D type) and $(D data).
     * $(D type) is checked and a $(D MessagePackException) is thrown if it does
     *  not match. The length of $(D data) is checked and a $(D MessagePackException)
     *  is thrown if the lengths do not match.  If $(D data) is null, a new slice
     *  is returned.
     */
    ref Unpacker unpackExt(ref byte type, ref ubyte[] data) return
    {
        import std.conv : text;

        canRead(Offset, 0);
        const header = read();

        uint length;
        uint rollbackLength = 0;
        if (header >= Format.EXT && header <= Format.EXT + 4)
        {
            // Fixed
            length = 2^^(header - Format.EXT);

        } else {
            // Dynamic length
            switch (header) with (Format)
            {
                case EXT8:
                    canRead(1);
                    length = read();
                    rollbackLength = 1;
                    break;
                case EXT16:
                    canRead(2);
                    length = load16To!ushort(read(2));
                    rollbackLength = 2;
                    break;
                case EXT32:
                    canRead(4);
                    length = load32To!uint(read(4));
                    rollbackLength = 4;
                    break;
                default:
                    rollback(0, "ext", cast(Format)header);
            }

        }

        canRead(1 + length);

        // Read and check the type
        byte type_ = read();
        rollbackLength += 1;
        if (type_ != type)
            rollback(rollbackLength, text("Cannot unpack EXT of type ", type_, " into type ", type));

        // Read and check data
        if (data is null)
            data = new ubyte[](length);
        else if (data.length != length) {
            rollback(rollbackLength, text("Length mismatch while unpacking EXT: ", data.length, " was given, actual length is ", length));
        }
        data[] = read(length);
        return this;
    }

    /**
     * Scans an entire buffer and converts each objects.
     *
     * This method is used for unpacking record-like objects.
     *
     * Example:
     * -----
     * // serialized data is "[1, 2][3, 4][5, 6][...".
     * auto unpacker = Unpacker(serializedData);
     * foreach (n, d; &unpacker.scan!(int, int))  // == "foreach (int n, int d; unpacker)"
     *     writeln(n, d); // 1st loop "1, 2", 2nd loop "3, 4"...
     * -----
     */
    int scan(Types...)(scope int delegate(ref Types) dg)
    {
        return opApply!(Types)(delegate int(ref Types objects) { return dg(objects); });
    }


    /// ditto
    int opApply(Types...)(scope int delegate(ref Types) dg)
    {
        int result;

        while (used_ - offset_) {
            auto length = beginArray();
            if (length != Types.length)
                rollback(calculateSize(length), "the number of deserialized objects is mismatched");

            Types objects;
            foreach (i, T; Types)
                unpack(objects[i]);

            result = dg(objects);
            if (result)
                return result;
        }

        return result;
    }


  private:
    /*
     * Deserializes nil object and assigns to $(D_PARAM value).
     *
     * Params:
     *  value = the reference of value to assign.
     *
     * Returns:
     *  self, i.e. for method chaining.
     *
     * Throws:
     *  UnpackException when doesn't read from buffer or precision loss occurs and
     *  MessagePackException when $(D_PARAM T) type doesn't match serialized type.
     */
    @safe
    ref Unpacker unpackNil(T)(ref T value)
    {
        canRead(Offset, 0);
        const header = read();

        if (header == Format.NIL)
            value = null;
        else
            rollback(0, "nil", cast(Format)header);

        return this;
    }


    /*
     * Next object is nil?
     *
     * Returns:
     *  true if next object is nil.
     */
    @safe
    bool checkNil()
    {
        canRead(Offset, 0);

        return buffer_[offset_] == Format.NIL;
    }


    /*
     * Calculates the format size of container length.
     */
    size_t calculateSize(bool rawType = false)(in size_t length)
    {
        static if (rawType)
            return length < 32 ? 0 : length < 65536 ? ushort.sizeof : uint.sizeof;
        else
            return length < 16 ? 0 : length < 65536 ? ushort.sizeof : uint.sizeof;
    }


    /*
     * Reading test.
     *
     * Params:
     *  size   = the size to read.
     *  offset = the offset to subtract when doesn't read from buffer.
     *
     * Throws:
     *  UnpackException when doesn't read from buffer.
     */
    @safe
    void canRead(in size_t size, in size_t offset = Offset)
    {
        if (used_ - offset_ < size) {
            if (offset)
                offset_ -= offset;

            throw new UnpackException("Insufficient buffer");
        }
    }


    /*
     * Reads value from buffer and advances offset.
     */
    @safe
    nothrow ubyte read()
    {
        return buffer_[offset_++];
    }


    /*
     * Reads value from buffer and advances offset.
     */
    @safe
    nothrow ubyte[] read(in size_t size)
    {
        auto result = buffer_[offset_..offset_ + size];

        offset_ += size;

        return result;
    }


    /*
     * Do rollback and throws exception.
     */
    @safe
    void rollback(in size_t size, in string reason)
    {
        offset_ -= size + Offset;
        onInvalidType(reason);
    }

    @safe
    void rollback(in size_t size, in string expected, in Format actual)
    {
        offset_ -= size + Offset;
        onInvalidType(expected, actual);
    }
}


private:


/*
 * A callback for type-mismatched error in deserialization process.
 */
@safe
pure void onInvalidType(in string reason)
{
    throw new MessagePackException("Attempt to unpack with non-compatible type: reason = " ~ reason);
}

@safe
pure void onInvalidType(in string expected, in Format actual)
{
    import std.conv: text;
    throw new MessagePackException(text("Attempt to unpack with non-compatible type: expected = ", expected, ", actual = ", actual));
}


unittest
{
    import msgpack.packer;

    { // unique
        mixin DefinePacker;

        Tuple!(bool, bool) result;
        Tuple!(bool, bool) test = tuple(true, false);

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(result);
        assert(test == result);
    }
    { // uint *
        mixin DefinePacker;

        Tuple!(ubyte, ushort, uint, ulong) result;
        Tuple!(ubyte, ushort, uint, ulong) test = tuple(cast(ubyte)ubyte.max, cast(ushort)ushort.max,
                                                        cast(uint)uint.max,   cast(ulong)ulong.max);

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(result);
        assert(test == result);
    }
    { // int *
        mixin DefinePacker;

        Tuple!(byte, short, int, long) result;
        Tuple!(byte, short, int, long) test = tuple(cast(byte)byte.min, cast(short)short.min,
                                                    cast(int)int.min,   cast(long)long.min);

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(result);
        assert(test == result);
    }
    { // floating point
        mixin DefinePacker;

        static if (real.sizeof == double.sizeof || !EnableReal)
        {
            Tuple!(float, double, double) result;
            Tuple!(float, double, double) test = tuple(cast(float)float.min_normal, cast(double)double.max, cast(real)double.min_normal);
        }
        else
        {
            Tuple!(float, double, real) result;
            Tuple!(float, double, real) test = tuple(cast(float)float.min_normal, cast(double)double.max, cast(real)real.min_normal);
        }

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(result);
        assert(test == result);
    }
    { // pointer
        mixin DefinePacker;

        Tuple!(ulong, long, double) origin;
        Tuple!(ulong, long, double) values = tuple(ulong.max, long.min, double.min_normal);
        Tuple!(ulong*, long*, double*) result = tuple(&origin.field[0], &origin.field[1], &origin.field[2]);
        Tuple!(ulong*, long*, double*) test = tuple(&values.field[0], &values.field[1], &values.field[2]);

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(result);
        foreach (i, v; test.field)
            assert(*v == *result.field[i]);
        assert(origin == values);
    }
    { // enum
        enum   : float { D = 0.5 }
        enum E : ulong { U = 100 }

        mixin DefinePacker;

        float f = D,   resultF;
        E     e = E.U, resultE;

        packer.pack(D, e);

        auto unpacker = Unpacker(packer.stream.data);

        unpacker.unpack(resultF, resultE);
        assert(f == resultF);
        assert(e == resultE);
    }
    { // container
        mixin DefinePacker;

        Tuple!(ulong[], double[uint], string, bool[2], char[2]) test
            = tuple([1UL, 2], [3U:4.0, 5:6.0, 7:8.0], "MessagePack is nice!", [true, false], "D!");

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);
        Tuple!(ulong[], double[uint], string, bool[2], char[2]) result;

        unpacker.unpack(result);
        assert(test == result);
    }
    { // ext

        // Try a variety of lengths, making sure to hit all the fixexts
        foreach (L; TypeTuple!(1, 2, 3, 4, 5, 8, 9, 16, 32, 512, 2^^16))
        {
            mixin DefinePacker;

            ubyte[] data = new ubyte[](L);
            ExtValue ext = ExtValue(7, data);
            packer.pack(ext);

            auto unpacker1 = Unpacker(packer.stream.data);
            ExtValue witness;

            unpacker1.unpack(witness);
            assert(ext == witness);

            // And try unpackExt
            auto unpacker2 = Unpacker(packer.stream.data);
            byte type = 1;
            ubyte[] deserializedData = new ubyte[](7);

            // This should be a type mismatch (1 != 7)
            assertThrown!MessagePackException(
                unpacker2.unpackExt(type, deserializedData));
            type = 7;

            // A data size mismatch
            assertThrown!MessagePackException(
                unpacker2.unpackExt(type, deserializedData));
            deserializedData = new ubyte[](L);

            // And this should succeed
            unpacker2.unpackExt(type, deserializedData);
            assert(deserializedData == data);
        }
    }
    { // user defined
        {
            static struct S
            {
                uint num;

                void toMsgpack(P)(ref P p) const { p.packArray(num); }
                void fromMsgpack(ref Unpacker u)
                {
                    assert(u.beginArray() == 1);
                    u.unpack(num);
                }
            }

            mixin DefinePacker; S result, test = S(uint.max);

            packer.pack(test);

            auto unpacker = Unpacker(packer.stream.data);
            unpacker.unpack(result);

            assert(test.num == result.num);
        }
        {
            static class C
            {
                uint num;

                this(uint n) { num = n; }

                void toMsgpack(P)(ref P p) const { p.packArray(num - 1); }
                void fromMsgpack(ref Unpacker u)
                {
                    assert(u.beginArray() == 1);
                    u.unpack(num);
                }
            }

            mixin DefinePacker; C result, test = new C(ushort.max);

            packer.pack(test);

            auto unpacker = Unpacker(packer.stream.data);
            unpacker.unpack(result, ushort.max);

            assert(test.num == result.num + 1);
        }
    }
    { // simple struct and class
        {
            static struct Simple
            {
                uint num;
                @nonPacked string str;
            }

            static struct Simple2
            {
                @nonPacked string str;
                uint num;
            }

            foreach (Type; TypeTuple!(Simple, Simple2)) {
                mixin DefinePacker;
                Type result, test;
                test.num = uint.max;
                test.str = "ignored";

                packer.pack(test);
                auto unpacker = Unpacker(packer.stream.data);
                unpacker.unpack(result);

                assert(test.num == result.num);
                assert(test.str != result.str);
            }
        }

        {
            static struct SimpleProxy1
            {
                import std.conv;
                static void serialize(ref Packer p, ref string val) { p.pack(to!uint(val)); }
                static void deserialize(ref Unpacker u, ref string val) { uint tmp; u.unpack(tmp); val = to!string(tmp); }
            }
            static struct SimpleWithProxied1
            {
                @serializedAs!SimpleProxy1 string data;
                enum string defaultValue = "10";
            }

            // https://github.com/msgpack/msgpack-d/issues/83
            static struct SimpleProxy2
            {
                import std.datetime;
                static void serialize(ref Packer p, ref SysTime val) { p.pack(val.toISOExtString()); }
                static void deserialize(ref Unpacker u, ref SysTime val) { string tmp; u.unpack(tmp); val = SysTime.fromISOExtString(tmp); }
            }
            static struct SimpleWithProxied2
            {
                import std.datetime;
                @serializedAs!SimpleProxy2 SysTime data;
                static SysTime defaultValue() @property { return SysTime(DateTime(2019,1,1,0,0,0)); }
            }
            
            foreach (Type; TypeTuple!(SimpleWithProxied1, SimpleWithProxied2)) {
                mixin DefinePacker;
                Type result, test;
                test.data = Type.defaultValue;
                
                packer.pack(test);
                auto unpacker = Unpacker(packer.stream.data);
                unpacker.unpack(result);
                assert(test.data == result.data);
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
            @nonPacked string str;
        }

        static class SimpleC2 : SimpleB
        {
            @nonPacked string str;
            uint num = uint.max;
        }

        static class SimpleD
        {
            static struct Proxy
            {
                import std.conv;
                static void serialize(ref Packer p, ref bool val)      { p.pack(to!string(val)); }
                static void serialize(ref Packer p, ref uint val)      { p.pack(to!string(val)); }
                static void serialize(ref Packer p, ref ubyte val)     { p.pack(to!string(val)); }
                static void deserialize(ref Unpacker u, ref bool val)  { string tmp; u.unpack(tmp); val = to!bool(tmp); }
                static void deserialize(ref Unpacker u, ref uint val)  { string tmp; u.unpack(tmp); val = to!uint(tmp); }
                static void deserialize(ref Unpacker u, ref ubyte val) { string tmp; u.unpack(tmp); val = to!ubyte(tmp); }
            }
            @serializedAs!Proxy bool flag = true;
            @serializedAs!Proxy ubyte type = 100;
            @serializedAs!Proxy uint num = uint.max;
            @nonPacked string str;
        }

        { // from derived class
            foreach (Type; TypeTuple!(SimpleC, SimpleC2, SimpleD)) {
                mixin DefinePacker;
                Type result, test = new Type();
                test.flag = false;
                test.type = 99;
                test.num  = uint.max / 2;
                test.str  = "ignored";

                packer.pack(test);
                auto unpacker = Unpacker(packer.stream.data);
                unpacker.unpack(result);

                assert(test.flag == result.flag);
                assert(test.type == result.type);
                assert(test.num  == result.num);
                assert(test.str  != result.str);
            }
        }
        { // from base class
            mixin DefinePacker; SimpleC test = new SimpleC();

            packer.pack(test);

            SimpleB result = new SimpleC();
            auto unpacker  = Unpacker(packer.stream.data);

            try {
                unpacker.unpack(result);
                assert(false);
            } catch (Exception e) { }
        }
        { // https://github.com/msgpack/msgpack-d/issues/16
            mixin DefinePacker;

            static class Issue16
            {
                int i;
                this(int i) { this.i = i; }
            }

            Issue16 c1 = new Issue16(10);

            // change behaviour to accept null with new object without constructor
            Issue16 c2 = null;
            packer.pack(c1);
            auto unpacker1 = Unpacker(packer.stream.data);
            unpacker1.unpack(c2);
            //unpack(pack(c1), c2);
            assert(c2.i == c1.i);

            Issue16 c3 = new Issue16(20);
            packer.stream.clear();
            packer.pack(c1);
            auto unpacker2 = Unpacker(packer.stream.data);
            unpacker2.unpack(c3);
            //unpack(pack(c1), c3);
            assert(c3.i == c1.i);
        }
    }
    { // variadic
        mixin DefinePacker;

        Tuple!(uint, long, double) test = tuple(uint.max, long.min, double.max);

        packer.pack(test);

        auto unpacker = Unpacker(packer.stream.data);

        uint u; long l; double d;

        unpacker.unpackArray(u, l, d);
        assert(test == tuple(u, l, d));
    }
    { // scan / opApply
        ubyte[] data;
        mixin DefinePacker;

        foreach (i; 0..2)
            packer.pack(tuple(1, 0.5, "Hi!"));

        foreach (n, d, s; &Unpacker(packer.stream.data).scan!(int, double, string)) {
            assert(n == 1);
            assert(d == 0.5);
            assert(s == "Hi!");
        }
    }
}
