module msgpack.value;

import msgpack.common;
import msgpack.attribute;
import msgpack.exception;

import std.json;
import std.container : Array;
import std.traits;
import std.typecons : Tuple, isTuple;


/**
 * $(D Value) is a $(D MessagePack) value representation
 *
 * Example:
 * -----
 * auto unpacker = StreamingUnpacker(pack(1, 0.1L) ~ pack(true) ~ pack("foobarbaz"));
 *
 * foreach (unpacked; unpacker) {
 *     if (unpacked.type == Value.Type.array) {
 *         foreach (obj; unpacked) {
 *             switch (obj.type) {
 *             case Value.Type.unsigned: writeln(obj.as!(uint)); break;
 *             case Value.Type.floating:            writeln(obj.as!(real)); break;
 *             defalut:
 *                 throw new Exception("Unknown type");
 *             }
 *         }
 *     } else {
 *         if (unpacked.type == Value.Type.boolean)
 *             writeln(unpacked.as!(bool));
 *         else
 *             writeln("Message: ", unpacked.as!(string));
 *     }
 * }
 * -----
 */
struct Value
{
    /**
     * $(D MessagePack) value type
     */
    static enum Type
    {
        nil,       /// nil(null in D)
        boolean,   /// true, false
        unsigned,  /// positive fixnum, uint 8, uint 16, uint 32, uint 64
        signed,    /// negative fixnum, int 8, int 16, int 32, int 64
        floating,  /// float, double, real
        array,     /// fix array, array 16, array 32
        map,       /// fix map, map 16, map 32
        raw,       /// fix raw, raw 16, raw 32
        ext        /// fix ext, ext8, ext16, ext32
    }


    /**
     * msgpack value representation
     */
    static union Via
    {
        bool         boolean;   /// corresponding to Type.boolean
        ulong        uinteger;  /// corresponding to Type.unsigned
        long         integer;   /// corresponding to Type.signed
        real         floating;  /// corresponding to Type.floating
        Value[]      array;     /// corresponding to Type.array
        Value[Value] map;       /// corresponding to Type.map
        ubyte[]      raw;       /// corresponding to Type.raw
        ExtValue     ext;       /// corresponding to Type.ext
    }


    Type type;  /// represents value type
    Via  via;   /// represents real value


    /**
     * Constructs a $(D Value) with arguments.
     *
     * Params:
     *  value = the real content.
     *  type  = the type of value.
     */
    @safe
    this(Type type)
    {
        this.type = type;
    }

    @safe
    this(typeof(null))
    {
        this(Type.nil);
    }

    /// ditto
    @trusted
    this(bool value, Type type = Type.boolean)
    {
        this(type);
        via.boolean = value;
    }


    /// ditto
    @trusted
    this(ulong value, Type type = Type.unsigned)
    {
        this(type);
        via.uinteger = value;
    }


    /// ditto
    @trusted
    this(long value, Type type = Type.signed)
    {
        this(type);
        via.integer = value;
    }


    /// ditto
    @trusted
    this(real value, Type type = Type.floating)
    {
        this(type);
        via.floating = value;
    }


    /// ditto
    @trusted
    this(Value[] value, Type type = Type.array)
    {
        this(type);
        via.array = value;
    }


    /// ditto
    @trusted
    this(Value[Value] value, Type type = Type.map)
    {
        this(type);
        via.map = value;
    }


    /// ditto
    @trusted
    this(ubyte[] value, Type type = Type.raw)
    {
        this(type);
        via.raw = value;
    }

    /// This is unsafe overload because using cast internally.
    @trusted
    this(string value, Type type = Type.raw)
    {
        this(type);
        via.raw = cast(ubyte[])value;
    }

    /**
     * Constructs a $(D Value) with arguments.
     *
     * Params:
     *  value = the real content.
     *  type  = the type of value.
     */
    @trusted
    this(ExtValue value, Type type = Type.ext)
    {
        this(type);
        via.ext = value;
    }

    /**
     * Converts value to $(D_PARAM T) type.
     *
     * Returns:
     *  converted value.
     *
     * Throws:
     *  MessagePackException if type is mismatched.
     *
     * NOTE:
     *  Current implementation uses cast.
     */
    @property @trusted
    T as(T)() if (is(Unqual!T == bool))
    {
        if (type != Type.boolean)
            onCastError();

        return via.boolean;
    }


    /// ditto
    @property @trusted
    T as(T)() if (isIntegral!T && !is(Unqual!T == enum))
    {
        if (type == Type.unsigned)
            return cast(T)via.uinteger;

        if (type == Type.signed)
            return cast(T)via.integer;

        onCastError();

        assert(false);
    }


    /// ditto
    @property @trusted
    T as(T)() if (isFloatingPoint!T && !is(Unqual!T == enum))
    {
        if (type != Type.floating)
            onCastError();

        return cast(T)via.floating;
    }


    /// ditto
    @property @trusted
    T as(T)() if (is(Unqual!T == enum))
    {
        return cast(T)as!(OriginalType!T);
    }


    /// ditto
    @property @trusted
    T as(T)() if (is(Unqual!T == ExtValue))
    {
        if (type != Type.ext)
            onCastError();

        return cast(T)via.ext;
    }


    /// ditto
    @property @trusted
    T as(T)() if ((isArray!T ||
                   isInstanceOf!(Array, T)) &&
                  !is(Unqual!T == enum))
    {
        alias typeof(T.init[0]) V;

        if (type == Type.nil) {
            static if (isDynamicArray!T) {
                return null;
            } else {
                return T.init;
            }
        }

        static if (isByte!V || isSomeChar!V) {
            if (type != Type.raw)
                onCastError();

            static if (isDynamicArray!T) {
                return cast(T)via.raw;
            } else {
                if (via.raw.length != T.length)
                    onCastError();

                return cast(T)(via.raw[0 .. T.length]);
            }
        } else {
            if (type != Type.array)
                onCastError();

            V[] array;

            foreach (elem; via.array)
                array ~= elem.as!(V);

            return array;
        }
    }


    /// ditto
    @property @trusted
    T as(T)() if (isAssociativeArray!T)
    {
        alias typeof(T.init.keys[0])   K;
        alias typeof(T.init.values[0]) V;

        if (type == Type.nil)
            return null;

        if (type != Type.map)
            onCastError();

        V[K] map;

        foreach (key, value; via.map)
            map[key.as!(K)] = value.as!(V);

        return map;
    }


    /**
     * Converts to $(D_PARAM T) type.
     *
     * Calling $(D fromMsgpack) if $(D_KEYWORD class) and $(D_KEYWORD struct) implement $(D fromMsgpack) method. $(D fromMsgpack) signature is:
     * -----
     * void fromMsgpack(Value value)
     * -----
     * This method assigns converted values to all members of T object if $(D_KEYWORD class) and $(D_KEYWORD struct) don't implement $(D fromMsgpack).
     *
     * Params:
     *  args = arguments to class constructor(class only).
     *
     * Returns:
     *  converted value.
     */
    @property @trusted
    T as(T, Args...)(Args args) if (is(Unqual!T == class))
    {
        if (type == Type.nil)
            return null;

        T object = new T(args);

        static if (hasMember!(T, "fromMsgpack"))
        {
            static if (__traits(compiles, { object.fromMsgpack(this); })) {
                object.fromMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'fromMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else {
            alias SerializingClasses!(T) Classes;

            if (via.array.length != SerializingMemberNumbers!(Classes))
                throw new MessagePackException("The number of deserialized object member is mismatched");

            size_t offset;
            foreach (Class; Classes) {
                Class obj = cast(Class)object;
                foreach (i, member; obj.tupleof) {
                    static if (isPackedField!(Class.tupleof[i]))
                        obj.tupleof[i] = via.array[offset++].as!(typeof(member));
                }
            }
        }

        return object;
    }


    /// ditto
    @property @trusted
    T as(T)() if (is(Unqual!T == struct) && !is(Unqual!T == ExtValue))
    {
        T obj;

        static if (hasMember!(T, "fromMsgpack"))
        {
            static if (__traits(compiles, { obj.fromMsgpack(this); })) {
                obj.fromMsgpack(this);
            } else {
                static assert(0, "Failed to invoke 'fromMsgpack' on type '" ~ Unqual!T.stringof ~ "'");
            }
        } else {
            static if (isTuple!T) {
                if (via.array.length != T.Types.length)
                    throw new MessagePackException("The number of deserialized Tuple element is mismatched");

                foreach (i, Type; T.Types)
                    obj.field[i] = via.array[i].as!(Type);
            } else {  // simple struct
                if (via.array.length != SerializingMemberNumbers!T)
                    throw new MessagePackException("The number of deserialized struct member is mismatched");

                size_t offset;
                foreach (i, member; obj.tupleof) {
                    static if (isPackedField!(T.tupleof[i]))
                        obj.tupleof[i] = via.array[offset++].as!(typeof(member));
                }
            }
        }

        return obj;
    }


    /**
     * Special method called by $(D Packer).
     *
     * Params:
     *  packer = a MessagePack serializer.
     */
    void toMsgpack(Packer)(ref Packer packer) const
    {
        final switch (type) {
        case Type.nil:
            packer.pack(null);
            break;
        case Type.boolean:
            packer.pack(via.boolean);
            break;
        case Type.unsigned:
            packer.pack(via.uinteger);
            break;
        case Type.signed:
            packer.pack(via.integer);
            break;
        case Type.floating:
            packer.pack(via.floating);
            break;
        case Type.raw:
            packer.pack(via.raw);
            break;
        case Type.ext:
            packer.packExt(via.ext.type, via.ext.data);
            break;
        case Type.array:
            packer.beginArray(via.array.length);
            foreach (elem; via.array)
                elem.toMsgpack(packer);
            break;
        case Type.map:
            packer.beginMap(via.map.length);
            foreach (key, value; via.map) {
                key.toMsgpack(packer);
                value.toMsgpack(packer);
            }
            break;
        }
    }


    /**
     * Comparison for equality. @trusted for union.
     */
    @trusted
    bool opEquals(Tdummy = void)(ref const Value other) const
    {
        if (type != other.type)
            return false;

        final switch (other.type) {
        case Type.nil:      return true;
        case Type.boolean:  return opEquals(other.via.boolean);
        case Type.unsigned: return opEquals(other.via.uinteger);
        case Type.signed:   return opEquals(other.via.integer);
        case Type.floating: return opEquals(other.via.floating);
        case Type.raw:      return opEquals(other.via.raw);
        case Type.ext:      return opEquals(other.via.ext);
        case Type.array:    return opEquals(other.via.array);
        case Type.map:      return opEquals(other.via.map);
        }
    }


    /// ditto
    @trusted
    bool opEquals(T : bool)(in T other) const
    {
        if (type != Type.boolean)
            return false;

        return via.boolean == other;
    }


    /// ditto
    @trusted
    bool opEquals(T : ulong)(in T other) const
    {
        static if (__traits(isUnsigned, T)) {
            if (type != Type.unsigned)
                return false;

            return via.uinteger == other;
        } else {
            if (type != Type.signed)
                return false;

            return via.integer == other;
        }
    }


    /// ditto
    @trusted
    bool opEquals(T : real)(in T other) const
    {
        if (type != Type.floating)
            return false;

        return via.floating == other;
    }


    /// ditto
    @trusted
    bool opEquals(T : const Value[])(in T other) const
    {
        if (type != Type.array)
            return false;

        return via.array == other;
    }


    /// ditto
    @trusted
    bool opEquals(T : const Value[Value])(in T other) const
    {
        if (type != Type.map)
            return false;

        // This comparison is instead of default comparison because 'via.map == other' raises "Access Violation".
        foreach (key, value; via.map) {
            if (key in other) {
                if (other[key] != value)
                    return false;
            } else {
                return false;
            }
        }

        return true;
    }


    /// ditto
    @trusted
    bool opEquals(T : const(ubyte)[])(in T other) const
    {
        if (type != Type.raw)
            return false;

        return via.raw == other;
    }


    /// ditto
    @trusted
    bool opEquals(T : string)(in T other) const
    {
        if (type != Type.raw)
            return false;

        return via.raw == cast(ubyte[])other;
    }


    //
    @trusted
    bool opEquals(T : ExtValue)(in T other) const
    {
        if (type != Type.ext)
            return false;

        return via.ext.type == other.type && via.ext.data == other.data;
    }


    @trusted
    hash_t toHash() const nothrow
    {
        static hash_t getHash(T)(T* v) @safe nothrow
        {
            return typeid(T).getHash(v);
        }

        final switch (type) {
        case Type.nil:      return 0;
        case Type.boolean:  return getHash(&via.boolean);
        case Type.unsigned: return getHash(&via.uinteger);
        case Type.signed:   return getHash(&via.integer);
        case Type.floating: return getHash(&via.floating);
        case Type.raw:      return getHash(&via.raw);
        case Type.ext:      return getHash(&via.ext);
        case Type.array:
            hash_t ret;
            foreach (elem; via.array)
                ret ^= elem.toHash();
            return ret;
        case Type.map:
            try {
                hash_t ret;
                foreach (key, value; via.map) {
                    ret ^= key.toHash();
                    ret ^= value.toHash();
                }
                return ret;
            } catch(Throwable) assert(0);
        }
    }
}


unittest
{
    import std.array;

    // nil
    Value value = Value(null);
    Value other = Value();

    assert(value      == other);
    assert(value.type == Value.Type.nil);

    // boolean
    value = Value(true);
    other = Value(false);

    assert(value           != other);
    assert(value.type      == Value.Type.boolean);
    assert(value.as!(bool) == true);
    assert(other           == false);

    try {
        auto b = value.as!(uint);
        assert(false);
    } catch (MessagePackException e) { }

    // unsigned integer
    value = Value(10UL);
    other = Value(10UL);

    assert(value           == other);
    assert(value.type      == Value.Type.unsigned);
    assert(value.as!(uint) == 10);
    assert(other           == 10UL);

    // signed integer
    value = Value(-20L);
    other = Value(-10L);

    assert(value          != other);
    assert(value.type     == Value.Type.signed);
    assert(value.as!(int) == -20);
    assert(other          == -10L);

    // enum
    enum E : int { F = -20 }

    E e = value.as!(E);
    assert(e == E.F);

    // floating point
    value = Value(0.1e-10L);
    other = Value(0.1e-20L);

    assert(value           != other);
    assert(value.type      == Value.Type.floating);
    assert(value.as!(real) == 0.1e-10L);
    assert(other           == 0.1e-20L);

    // raw
    value = Value(cast(ubyte[])[72, 105, 33]);
    other = Value(cast(ubyte[])[72, 105, 33]);

    assert(value               == other);
    assert(value.type          == Value.Type.raw);
    assert(value.as!(string)   == "Hi!");
    assert(value.as!(ubyte[3]) == [72, 105, 33]);
    assert(other               == cast(ubyte[])[72, 105, 33]);

    // raw with string
    value = Value("hello");
    other = Value("hello");

    assert(value             == other);
    assert(value.type        == Value.Type.raw);
    assert(value.as!(string) == "hello");

    // enum : string
    enum EStr : string { elem = "hello" }

    assert(value.as!(EStr) == EStr.elem);

    // ext
    auto ext = ExtValue(7, [1,2,3]);
    value = Value(ExtValue(7, [1,2,3]));
    assert(value.as!ExtValue == ext);

    // array
    auto t = Value(cast(ubyte[])[72, 105, 33]);
    value = Value([t]);
    other = Value([t]);

    assert(value               == other);
    assert(value.type          == Value.Type.array);
    assert(value.as!(string[]) == ["Hi!"]);
    assert(other               == [t]);

    // map
    value = Value([Value(1L):Value(2L)]);
    other = Value([Value(1L):Value(1L)]);

    assert(value               != other);
    assert(value.type          == Value.Type.map);
    assert(value.as!(int[int]) == [1:2]);
    assert(other               == [Value(1L):Value(1L)]);

    value = Value(10UL);

    // struct
    static struct S
    {
        ulong num;

        void fromMsgpack(Value value) { num = value.via.uinteger; }
    }

    S s = value.as!(S);
    assert(s.num == 10);

    value = Value([Value(0.5f), Value(cast(ubyte[])[72, 105, 33])]);

    // struct
    static struct Simple
    {
        @nonPacked int era;
        double num;
        string msg;
    }

    Simple simple = value.as!(Simple);
    assert(simple.era == int.init);
    assert(simple.num == 0.5f);
    assert(simple.msg == "Hi!");

    value = Value(10UL);

    // class
    static class C
    {
        ulong num;

        void fromMsgpack(Value value) { num = value.via.uinteger; }
    }

    C c = value.as!(C);
    assert(c.num == 10);

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
        @nonPacked string str;
        uint num = uint.max;
    }

    value = Value([Value(false), Value(99UL), Value(cast(ulong)(uint.max / 2u))]);

    SimpleC sc = value.as!(SimpleC);
    assert(sc.flag == false);
    assert(sc.type == 99);
    assert(sc.num  == uint.max / 2);
    assert(sc.str.empty);

    // std.typecons.Tuple
    value = Value([Value(true), Value(1UL), Value(cast(ubyte[])"Hi!")]);

    auto tuple = value.as!(Tuple!(bool, uint, string));
    assert(tuple.field[0] == true);
    assert(tuple.field[1] == 1u);
    assert(tuple.field[2] == "Hi!");

    /*
     * non-MessagePackable object is stopped by static assert
     * static struct NonMessagePackable {}
     * auto nonMessagePackable = value.as!(NonMessagePackable);
     */
}


/**
 * Converts $(D Value) to $(D JSONValue).
 *
 * Params:
 *  val = $(D Value) to convert.
 *
 * Returns:
 *  a $(D JSONValue).
 */
@trusted
JSONValue toJSONValue(in Value val)
{
    final switch (val.type)
    {
        case Value.Type.nil:      return JSONValue(null);
        case Value.Type.boolean:  return JSONValue(val.via.boolean);
        case Value.Type.unsigned: return JSONValue(val.via.uinteger);
        case Value.Type.signed:   return JSONValue(val.via.integer);
        case Value.Type.floating: return JSONValue(val.via.floating);
        case Value.Type.raw:      return JSONValue(cast(string)(val.via.raw.idup));
        case Value.Type.ext:      throw new MessagePackException("Unable to convert ext to json");
        case Value.Type.array: {
            JSONValue[] vals;
            foreach (elem; val.via.array)
                vals ~= elem.toJSONValue();
            return JSONValue(vals);
        }
        case Value.Type.map: {
            JSONValue[string] vals;
            foreach (key, value; val.via.map) {
                if (key.type != Value.Type.raw)
                {
                    throw new MessagePackException("JSON-object key must be a raw type");
                }
                vals[key.as!string] = value.toJSONValue();
            }
            return JSONValue(vals);
        }
    }
}

/**
 * Converts $(D JSONValue) to $(D Value).
 *
 * Params:
 *  val = $(D JSONValue) to convert.
 *
 * Returns:
 *  a $(D Value).
 */
@trusted
Value fromJSONValue(in JSONValue val)
{
    final switch (val.type())
    {
        case JSONType.null_:      return Value(null);
        case JSONType.true_:      return Value(true);
        case JSONType.false_:     return Value(false);
        case JSONType.uinteger:  return Value(val.uinteger);
        case JSONType.integer:   return Value(val.integer);
        case JSONType.float_:     return Value(val.floating);
        case JSONType.string:    return Value(cast(ubyte[])(val.str));
        case JSONType.array: {
            Value[] vals;
            foreach (elem; val.array)
                vals ~= elem.fromJSONValue();
            return Value(vals);
        }
        case JSONType.object: {
            Value[Value] vals;
            foreach (key, value; val.object) {
                vals[Value(cast(ubyte[])key)] = value.fromJSONValue();
            }
            return Value(vals);
        }
    }
}

unittest
{
    import std.array : array;
    import std.algorithm : equal, map;
    import std.conv;
    import std.math : isClose;
    import std.range;
    import msgpack;

    // nil
    Value value = Value(null);

    assert(toJSONValue(value).type() == JSONType.null_);

    // boolean
    value = Value(true);
    auto other = Value(false);

    assert(toJSONValue(value).type() == JSONType.true_);
    assert(toJSONValue(other).type() == JSONType.false_);

    // unsigned integer
    value = Value(10UL);

    assert(value.toJSONValue().type == JSONType.uinteger);
    assert(value.toJSONValue().uinteger == value.as!uint);
    assert(value.toJSONValue().uinteger == 10UL);

    // signed integer
    value = Value(-20L);

    assert(value.toJSONValue().type == JSONType.integer);
    assert(value.toJSONValue().integer == value.as!int);

    // enum
    enum E : int { F = -20 }
    value = Value(cast(long)(E.F));

    assert(value.toJSONValue().type == JSONType.integer);
    assert(value.toJSONValue().integer == E.F);

    // floating point
    value = Value(0.1e-10L);
    other = Value(0.1e-20L);

    assert(value.toJSONValue().type == JSONType.float_);
    assert(other.toJSONValue().type == JSONType.float_);

    assert(isClose(value.toJSONValue().floating, 0.1e-10L));
    assert(isClose(other.toJSONValue().floating, 0.1e-20L));

    // raw
    long[] arr = [72, 105, 33];
    value = Value(to!(ubyte[])(arr));

    assert(value.toJSONValue().type == JSONType.string);
    assert(equal(value.toJSONValue().str, arr));

    // raw with string
    value = Value("hello");
    assert(value.toJSONValue().type == JSONType.string);
    assert(value.toJSONValue().str == "hello");

    // array
    auto t = Value(to!(ubyte[])(arr));
    value = Value([t]);
    other = Value(array(map!(a => Value(a))(arr)));

    assert(value.toJSONValue().type == JSONType.array);
    assert(value.toJSONValue().array.length == 1);
    assert(value.toJSONValue().array.front().type == JSONType.string);
    assert(equal(value.toJSONValue().array.front().str, arr));
    assert(other.toJSONValue().type == JSONType.array);
    assert(array(map!(a => a.integer)(other.toJSONValue().array)) == arr);

    // map
    value = Value([Value("key"):Value(2L)]);

    assert(value.toJSONValue().type == JSONType.object);
    assert("key" in value.toJSONValue().object);
    assert(value.toJSONValue().object["key"].type == JSONType.integer);
    assert(value.toJSONValue().object["key"].integer == 2L);

    // struct
    static struct Simple
    {
        @nonPacked int era;
        double num;
        string msg;
    }

    Simple simple;
    simple.era = 5;
    simple.num = 13.5;
    simple.msg = "helloworld";
    value = simple.pack().unpack().value;

    assert(value.toJSONValue().type == JSONType.array);
    assert(value.toJSONValue().array.length == 2);
    assert(value.toJSONValue().array[0].type == JSONType.float_);
    assert(isClose(value.toJSONValue().array[0].floating, simple.num));
    assert(value.toJSONValue().array[1].type == JSONType.string);
    assert(value.toJSONValue().array[1].str == simple.msg);

    // class
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
        @nonPacked string str;
        uint num = uint.max;
    }

    SimpleC sc = new SimpleC;
    value = sc.pack!true().unpack().value;

    assert(value.toJSONValue().type == JSONType.object);
    assert(value.toJSONValue().object.length == 3);
    assert("flag" in value.toJSONValue().object);
    assert(value.toJSONValue().object["flag"].type == (sc.flag ? JSONType.true_ : JSONType.false_));
    assert("type" in value.toJSONValue().object);
    assert(value.toJSONValue().object["type"].type == JSONType.uinteger);
    assert(value.toJSONValue().object["type"].uinteger == sc.type);
    assert("num" in value.toJSONValue().object);
    assert(value.toJSONValue().object["num"].type == JSONType.uinteger);
    assert(value.toJSONValue().object["num"].uinteger == sc.num);

    other = value.toJSONValue().fromJSONValue();
    assert(value == other);
}


private:


/**
 * A callback for type-mismatched error in cast conversion.
 */
@safe
pure void onCastError()
{
    throw new MessagePackException("Attempt to cast with another type");
}
