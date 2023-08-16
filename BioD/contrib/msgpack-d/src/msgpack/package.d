// Written in the D programming language.

/**
 * MessagePack serializer and deserializer implementation.
 *
 * MessagePack is a binary-based serialization specification.
 *
 * Example:
 * -----
 * auto data = tuple("MessagePack!", [1, 2], true);
 *
 * auto serialized = pack(data);
 *
 * // ...
 *
 * typeof(data) deserialized;
 *
 * unpack(serialized, deserialized);
 *
 * assert(data == deserialized);
 * -----
 *
 * See_Also:
 *  $(LINK2 http://msgpack.org/, The MessagePack Project)$(BR)
 *  $(LINK2 https://github.com/msgpack/msgpack/blob/master/spec.md, MessagePack data format)
 *
 * Copyright: Copyright Masahiro Nakagawa 2010-.
 * License:   <a href="http://www.boost.org/LICENSE_1_0.txt">Boost License 1.0</a>.
 * Authors:   Masahiro Nakagawa
 */
module msgpack;

public:

import msgpack.common;
import msgpack.attribute;
import msgpack.buffer;
import msgpack.exception;
import msgpack.packer;
import msgpack.unpacker;
import msgpack.streaming_unpacker;
import msgpack.register;
import msgpack.value;

version(Windows) {
    pragma(lib, "WS2_32");
}

@trusted:


/**
 * Serializes $(D_PARAM args).
 *
 * Assumes single object if the length of $(D_PARAM args) == 1,
 * otherwise array object.
 *
 * Params:
 *  args = the contents to serialize.
 *
 * Returns:
 *  a serialized data.
 */
ubyte[] pack(bool withFieldName = false, Args...)(in Args args)
{
    auto packer = Packer(withFieldName);

    static if (Args.length == 1)
        packer.pack(args[0]);
    else
        packer.packArray(args);

    return packer.stream.data;
}


/**
 * Deserializes $(D_PARAM buffer) using stream deserializer.
 *
 * Params:
 *  buffer = the buffer to deserialize.
 *
 * Returns:
 *  a $(D Unpacked) contains deserialized object.
 *
 * Throws:
 *  UnpackException if deserialization doesn't succeed.
 */
Unpacked unpack(in ubyte[] buffer)
{
    auto unpacker = StreamingUnpacker(buffer);

    if (!unpacker.execute())
        throw new UnpackException("Deserialization failure");

    return unpacker.unpacked;
}


/**
 * Deserializes $(D_PARAM buffer) using direct-conversion deserializer.
 *
 * Assumes single object if the length of $(D_PARAM args) == 1,
 * otherwise array object.
 *
 * Params:
 *  buffer = the buffer to deserialize.
 *  args   = the references of values to assign.
 */
void unpack(bool withFieldName = false, Args...)(in ubyte[] buffer, ref Args args)
{
    auto unpacker = Unpacker(buffer, buffer.length, withFieldName);

    static if (Args.length == 1)
        unpacker.unpack(args[0]);
    else
        unpacker.unpackArray(args);
}


/**
 * Return value version
 */
Type unpack(Type, bool withFieldName = false)(in ubyte[] buffer)
{
    auto unpacker = Unpacker(buffer, buffer.length, withFieldName);

    Type result;
    unpacker.unpack(result);
    return result;
}


unittest
{
    auto serialized = pack(false);

    assert(serialized[0] == Format.FALSE);

    auto deserialized = unpack(pack(1, true, "Foo"));

    assert(deserialized.type == Value.Type.array);
    assert(deserialized.via.array[0].type == Value.Type.unsigned);
    assert(deserialized.via.array[1].type == Value.Type.boolean);
    assert(deserialized.via.array[2].type == Value.Type.raw);
}


unittest
{
    import std.typecons;

    { // stream
        auto result = unpack(pack(false));

        assert(result.via.boolean == false);
    }
    { // direct conversion
        Tuple!(uint, string) result;
        Tuple!(uint, string) test = tuple(1, "Hi!");

        unpack(pack(test), result);
        assert(result == test);

        test.field[0] = 2;
        test.field[1] = "Hey!";
        unpack(pack(test.field[0], test.field[1]), result.field[0], result.field[1]);
        assert(result == test);
    }
    { // return value direct conversion
        Tuple!(uint, string) test = tuple(1, "Hi!");

        auto data = pack(test);
        assert(data.unpack!(Tuple!(uint, string)) == test);
    }
    { // serialize object as a Map
        static class C
        {
            int num;

            this(int num) { this.num = num; }
        }

        auto test = new C(10);
        auto result = new C(100);

        unpack!(true)(pack!(true)(test), result);
        assert(result.num == 10, "Unpacking with field names failed");
    }
}


unittest
{
    import std.typetuple;

    // unittest for https://github.com/msgpack/msgpack-d/issues/8
    foreach (Type; TypeTuple!(byte, short, int, long)) {
        foreach (i; [-33, -20, -1, 0, 1, 20, 33]) {
            Type a = cast(Type)i;
            Type b;
            unpack(pack(a), b);
            assert(a == b);
        }
    }
}


unittest
{
    import std.typetuple;

    // char types
    foreach (Type; TypeTuple!(char, wchar, dchar)) {
        foreach (i; [Type.init, Type.min, Type.max, cast(Type)'j']) {
            Type a = i;
            Type b;
            unpack(pack(a), b);
            assert(a == b);
        }
    }
}

unittest
{
    // ext type
    auto result = unpack(pack(ExtValue(7, [1,2,3,4])));
    assert(result == ExtValue(7, [1,2,3,4]));
}

unittest {
    import std.exception: assertThrown;

    struct Version {
        int major= -1;
        int minor = -1;
    }

    struct SubscriptionTopic {
        string[] topicComponents;
    }

    struct SubscriptionSender
    {
        string hostName;
        string biosName;
    }

    struct PubSubMessage {

        enum Type {
            publication,
            subscribe,
            unsubscribe,
        }

        Version version_;
        Type type;
        SubscriptionSender sender;
        SubscriptionTopic topic;
        string value;
    }

    ubyte[] bytes = [149, 146, 255, 255,   0, 146, 164, 104,
                     111, 115, 116, 164,  98, 105, 111, 115,
                     145, 221, 171, 105, 110, 116, 101, 114,
                     101, 115, 116, 105, 110, 103, 165, 116,
                     111, 112, 105,  99, 167, 112,  97, 121,
                     108, 111,  97, 100, 158, 142, 210,  31,
                     127,  81, 149, 125, 183, 108,  86,  17,
                     100,  35, 168];

    // should not throw OutOfMemoryError
    assertThrown!MessagePackException(unpack!PubSubMessage(bytes));
}


/**
 * Handy helper for creating MessagePackable object.
 *
 * toMsgpack / fromMsgpack are special methods for serialization / deserialization.
 * This template provides those methods to struct/class.
 *
 * Example:
 * -----
 * struct S
 * {
 *     int num; string str;
 *
 *     // http://d.puremagic.com/issues/show_bug.cgi?id = 1099
 *     mixin MessagePackable;  // all members
 *     // mixin MessagePackable!("num");  // num only
 * }
 * -----
 *
 * Defines those methods manually if you treat complex data-structure.
 */
mixin template MessagePackable(Members...)
{
    static if (Members.length == 0) {
        /**
         * Serializes members using $(D_PARAM packer).
         *
         * Params:
         *  packer = the serializer to pack.
         */
        void toMsgpack(Packer)(ref Packer packer, bool withFieldName = false) const
        {
            if (withFieldName) {
                packer.beginMap(this.tupleof.length);
                foreach (i, member; this.tupleof) {
                    packer.pack(getFieldName!(typeof(this), i));
                    packer.pack(member);
                }
            } else {
                packer.beginArray(this.tupleof.length);
                foreach (member; this.tupleof)
                    packer.pack(member);
            }
        }


        /**
         * Deserializes $(D MessagePack) object to members using Value.
         *
         * Params:
         *  value = the MessagePack value to unpack.
         *
         * Throws:
         *  MessagePackException if $(D_PARAM value) is not an Array type.
         */
        void fromMsgpack(Value value)
        {
            // enables if std.contracts.enforce is moved to object_.d
            // enforceEx!MessagePackException(value.type == Value.Type.array, "Value must be Array type");
            if (value.type != Value.Type.array)
                throw new MessagePackException("Value must be an Array type");
            if (value.via.array.length != this.tupleof.length)
                throw new MessagePackException("The size of deserialized value is mismatched");

            foreach (i, member; this.tupleof)
                this.tupleof[i] = value.via.array[i].as!(typeof(member));
        }


        /**
         * Deserializes $(D MessagePack) object to members using direct-conversion deserializer.
         *
         * Params:
         *  value = the reference to direct-conversion deserializer.
         *
         * Throws:
         *  MessagePackException if the size of deserialized value is mismatched.
         */
        void fromMsgpack(ref Unpacker unpacker)
        {
            auto length = unpacker.beginArray();
            if (length != this.tupleof.length)
                throw new MessagePackException("The size of deserialized value is mismatched");

            foreach (i, member; this.tupleof)
                unpacker.unpack(this.tupleof[i]);
        }
    } else {
        /**
         * Member selecting version of toMsgpack.
         */
        void toMsgpack(Packer)(ref Packer packer, bool withFieldName = false) const
        {
            if (withFieldName) {
                packer.beginMap(Members.length);
                foreach (member; Members) {
                    packer.pack(member);
                    packer.pack(mixin(member));
                }
            } else {
                packer.beginArray(Members.length);
                foreach (member; Members)
                    packer.pack(mixin(member));
            }
        }


        /**
         * Member selecting version of fromMsgpack for Value.
         */
        void fromMsgpack(Value value)
        {
            if (value.type != Value.Type.array)
                throw new MessagePackException("Value must be an Array type");
            if (value.via.array.length != Members.length)
                throw new MessagePackException("The size of deserialized value is mismatched");

            foreach (i, member; Members)
                mixin(member ~ "= value.via.array[i].as!(typeof(" ~ member ~ "));");
        }


        /**
         * Member selecting version of fromMsgpack for direct-converion deserializer.
         */
        void fromMsgpack(ref Unpacker unpacker)
        {
            auto length = unpacker.beginArray();
            if (length != Members.length)
                throw new MessagePackException("The size of deserialized value is mismatched");

            foreach (member; Members)
                unpacker.unpack(mixin(member));
        }
    }
}


unittest
{
    { // all members
        /*
         * Comment out because "src/msgpack.d(4048): Error: struct msgpack.__unittest16.S no size yet for forward reference" occurs
         */
        static struct S
        {
            uint num; string str;
            mixin MessagePackable;
        }

        mixin DefinePacker;

        S orig = S(10, "Hi!"); orig.toMsgpack(packer);

        { // stream
            auto unpacker = StreamingUnpacker(packer.stream.data); unpacker.execute();

            S result; result.fromMsgpack(unpacker.unpacked);

            assert(result.num == 10);
            assert(result.str == "Hi!");
        }
        { // direct conversion
            auto unpacker = Unpacker(packer.stream.data);

            S result; unpacker.unpack(result);

            assert(result.num == 10);
            assert(result.str == "Hi!");
        }
    }
    { // member select
        static class C
        {
            uint num; string str;

            this() {}
            this(uint n, string s) { num = n; str = s; }

            mixin MessagePackable!("num");
        }

        mixin DefinePacker;

        C orig = new C(10, "Hi!"); orig.toMsgpack(packer);

        { // stream
            auto unpacker = StreamingUnpacker(packer.stream.data); unpacker.execute();

            C result = new C; result.fromMsgpack(unpacker.unpacked);

            assert(result.num == 10);
        }
        { // direct conversion
            auto unpacker = Unpacker(packer.stream.data);

            C result; unpacker.unpack(result);

            assert(result.num == 10);
        }
    }
}


unittest
{
    import std.datetime: Clock, SysTime;
    import msgpack.packer, msgpack.unpacker;

    static struct SysTimePackProxy
    {
        static void serialize(ref Packer p, ref in SysTime tim)
        {
            p.pack(tim.toISOExtString());
        }

        static void deserialize(ref Unpacker u, ref SysTime tim)
        {
            string tmp;
            u.unpack(tmp);
            tim = SysTime.fromISOExtString(tmp);
        }
    }
    static struct LogData
    {
        string msg;
        string file;
        ulong  line;
        @serializedAs!SysTimePackProxy SysTime timestamp;

        this(string message, string file = __FILE__, ulong line = __LINE__)
        {
            this.msg = message;
            this.file = file;
            this.line = line;
            this.timestamp = Clock.currTime();
        }
    }

    /// Now we can serialize/deserialize LogData
    LogData[] logs;
    logs ~= LogData("MessagePack is nice!");
    auto data = pack(logs);
    LogData[] datas = unpack!(LogData[])(data);
    assert(datas[0].timestamp.toString() == datas[0].timestamp.toString());
}
