module msgpack.streaming_unpacker;

import msgpack.common;
import msgpack.attribute;
import msgpack.exception;
import msgpack.value;

import std.array;
import std.exception;
import std.range;
import std.stdio;
import std.traits;
import std.typecons;
import std.typetuple;
import std.container;


/**
 * $(D Unpacked) is a $(D Range) wrapper for stream deserialization result
 */
struct Unpacked
{
    import std.conv : text;

    Value value;  /// deserialized value

    alias value this;


    /**
     * Constructs a $(D Unpacked) with argument.
     *
     * Params:
     *  value = a deserialized value.
     */
    @safe
    this(ref Value value)
    {
        this.value = value;
    }


    /**
     * InputRange primitive operation that checks iteration state.
     *
     * Returns:
     *  true if there are no more elements to be iterated.
     */
    @property @trusted
    nothrow bool empty() const  // std.array.empty isn't nothrow function
    {
        return (value.type == Value.Type.array) && !value.via.array.length;
    }


    /**
     * Range primitive operation that returns the length of the range.
     *
     * Returns:
     *  the number of values.
     */
    @property @trusted
    size_t length()
    {
        debug enforce(value.type == Value.Type.array, "lenght is called with non array object. type = " ~ text(value.type));
        return value.via.array.length;
    }


    /**
     * InputRange primitive operation that returns the currently iterated element.
     *
     * Returns:
     *  the deserialized $(D Value).
     */
    @property @trusted
    ref Value front()
    {
        debug enforce(value.type == Value.Type.array, "front is called with non array object. type = " ~ text(value.type));
        return value.via.array.front;
    }


    /**
     * InputRange primitive operation that advances the range to its next element.
     */
    @trusted
    void popFront()
    {
        debug enforce(value.type == Value.Type.array, "popFront is called with non array object. type = " ~ text(value.type));
        value.via.array.popFront();
    }

    /**
     * RandomAccessRange primitive operation.
     *
     * Returns:
     *  the deserialized $(D Value) at $(D_PARAM n) position.
     */
    @trusted
    ref Value opIndex(size_t n)
    {
        debug enforce(value.type == Value.Type.array, "opIndex is called with non array object. type = " ~ text(value.type));
        return value.via.array[n];
    }

    /**
     * Returns a slice of the range.
     *
     * Paramas:
     *  from = the start point of slicing.
     *  to   = the end point of slicing.
     *
     * Returns:
     *  the slice of Values.
     */
    @trusted
    Value[] opSlice(size_t from, size_t to)
    {
        debug enforce(value.type == Value.Type.array, "opSlice is called with non array object. type = " ~ text(value.type));
        return value.via.array[from..to];
    }

    /**
     * Range primitive operation that returns the snapshot.
     *
     * Returns:
     *  the snapshot of this Value.
     */
    @property @safe
    Unpacked save()
    {
        return Unpacked(value);
    }
}


unittest
{
    static assert(isForwardRange!Unpacked);
    static assert(hasLength!Unpacked);
}


/**
 * This $(D StreamingUnpacker) is a $(D MessagePack) streaming deserializer
 *
 * This implementation enables you to load multiple objects from a stream(like network).
 *
 * Example:
 * -----
 * ...
 * auto unpacker = StreamingUnpacker(serializedData);
 * ...
 *
 * // appends new data to buffer if pre execute() call didn't finish deserialization.
 * unpacker.feed(newSerializedData);
 *
 * while (unpacker.execute()) {
 *     foreach (obj; unpacker.purge()) {
 *         // do stuff (obj is a Value)
 *     }
 * }
 *
 * if (unpacker.size)
 *     throw new Exception("Message is too large");
 * -----
 */
struct StreamingUnpacker
{
  private:
    /*
     * Context state of deserialization
     */
    enum State
    {
        HEADER = 0x00,

        BIN8 = 0x04,
        BIN16,
        BIN32,

        // Floating point, Unsigned, Signed interger (== header & 0x03)
        FLOAT = 0x0a,
        DOUBLE,
        UINT8,
        UINT16,
        UINT32,
        UINT64,
        INT8,
        INT16,
        INT32,
        INT64,

        // Container (== header & 0x01)
        STR8 = 0x19,
        RAW16 = 0x1a,
        RAW32,
        ARRAY16,
        ARRAY36,
        MAP16,
        MAP32,
        RAW,

        // EXT family
        EXT8,
        EXT16,
        EXT32,
        EXT_DATA,

        // D-specific type
        REAL
    }


    /*
     * Element type of container
     */
    enum ContainerElement
    {
        ARRAY_ITEM,
        MAP_KEY,
        MAP_VALUE
    }


    /*
     * Internal stack context
     */
    static struct Context
    {
        static struct Container
        {
            ContainerElement type;   // value container type
            Value            value;  // current value
            Value            key;    // for map value
            size_t           count;  // container length
        }

        State       state;  // current state of deserialization
        size_t      trail;  // current deserializing size
        size_t      top;    // current index of stack
        Container[] stack;  // storing values
    }

    Context context_;  // stack environment for streaming deserialization

    mixin InternalBuffer;


  public:
    /**
     * Constructs a $(D StreamingUnpacker).
     *
     * Params:
     *  target     = byte buffer to deserialize
     *  bufferSize = size limit of buffer size
     */
    @safe
    this(in ubyte[] target, in size_t bufferSize = 8192)
    {
        initializeBuffer(target, bufferSize);
        initializeContext();
    }


    /**
     * Forwards to deserialized object.
     *
     * Returns:
     *  the $(D Unpacked) object contains deserialized value.
     */
    @property @safe
    Unpacked unpacked()
    {
        return Unpacked(context_.stack[0].value);
    }


    /**
     * Clears some states for next deserialization.
     */
    @safe
    nothrow void clear()
    {
        initializeContext();

        parsed_ = 0;
    }


    /**
     * Convenient method for unpacking and clearing states.
     *
     * Example:
     * -----
     * foreach (obj; unpacker.purge()) {
     *     // do stuff
     * }
     * -----
     * is equivalent to
     * -----
     * foreach (obj; unpacker.unpacked) {
     *     // do stuff
     * }
     * unpacker.clear();
     * -----
     *
     * Returns:
     *  the $(D Unpacked) object contains deserialized value.
     */
    @safe
    Unpacked purge()
    {
        auto result = Unpacked(context_.stack[0].value);

        clear();

        return result;
    }


    /**
     * Executes deserialization.
     *
     * Returns:
     *  true if deserialization has been completed, otherwise false.
     *
     * Throws:
     *  $(D UnpackException) when parse error occurs.
     */
    bool execute()
    {
        /*
         * Current implementation is very dirty(goto! goto!! goto!!!).
         * This Complexity for performance(avoid function call).
         */

        bool   ret;
        size_t cur = offset_;
        Value obj;

        // restores before state
        auto state =  context_.state;
        auto trail =  context_.trail;
        auto top   =  context_.top;
        auto stack = &context_.stack;

        /*
         * Helper for container deserialization
         */
        void startContainer(string Type)(ContainerElement type, size_t length)
        {
            mixin("callback" ~ Type ~ "((*stack)[top].value, length);");

            (*stack)[top].type  = type;
            (*stack)[top].count = length;
            (*stack).length     = ++top + 1;
        }

        // non-deserialized data is nothing
        if (used_ - offset_ == 0)
            goto Labort;

        do {
          Lstart:
            if (state == State.HEADER) {
                const header = buffer_[cur];

                if (0x00 <= header && header <= 0x7f) {         // positive
                    callbackUInt(obj, header);
                    goto Lpush;
                } else if (0xe0 <= header && header <= 0xff) {  // negative
                    callbackInt(obj, cast(byte)header);
                    goto Lpush;
                } else if (0xa0 <= header && header <= 0xbf) {  // fix raw
                    trail = header & 0x1f;
                    state = State.RAW;
                    cur++;
                    continue;
                } else if (0xd4 <= header && header <= 0xd8) {  // fix ext
                    trail = 2 ^^ (header - 0xd4) + 1;
                    state = State.EXT_DATA;
                    cur++;
                    continue;
                } else if (0x90 <= header && header <= 0x9f) {  // fix array
                    size_t length = header & 0x0f;
                    if (length == 0) {
                        callbackArray(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Array"(ContainerElement.ARRAY_ITEM, length);
                        cur++;
                        continue;
                    }
                } else if (0x80 <= header && header <= 0x8f) {  // fix map
                    size_t length = header & 0x0f;
                    if (length == 0) {
                        callbackMap(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Map"(ContainerElement.MAP_KEY, length);
                        cur++;
                        continue;
                    }
                } else {
                    switch (header) {
                    case Format.UINT8, Format.UINT16, Format.UINT32, Format.UINT64,
                         Format.INT8, Format.INT16, Format.INT32, Format.INT64,
                         Format.FLOAT, Format.DOUBLE:
                        trail = 1 << (header & 0x03); // computes object size
                        state = cast(State)(header & 0x1f);
                        break;
                    case Format.REAL:
                        trail = RealSize;
                        state = State.REAL;
                        break;
                    case Format.ARRAY16, Format.ARRAY32,
                         Format.MAP16, Format.MAP32:
                        trail = 2 << (header & 0x01);  // computes container size
                        state = cast(State)(header & 0x1f);
                        break;
                    // raw will become str format in new spec
                    case Format.STR8:
                    case Format.RAW16: // will be STR16
                    case Format.RAW32: // will be STR32
                        trail = 1 << ((header & 0x03) - 1);  // computes container size
                        state = cast(State)(header & 0x1f);
                        break;
                    case Format.BIN8, Format.BIN16, Format.BIN32:
                        trail = 1 << (header & 0x03);  // computes container size
                        state = cast(State)(header & 0x1f);
                        break;
                    case Format.EXT8:
                        trail = 1;
                        state = State.EXT8;
                        break;
                    case Format.EXT16:
                        trail = 2;
                        state = State.EXT16;
                        break;
                    case Format.EXT32:
                        trail = 4;
                        state = State.EXT32;
                        break;
                    case Format.NIL:
                        callbackNil(obj);
                        goto Lpush;
                    case Format.TRUE:
                        callbackBool(obj, true);
                        goto Lpush;
                    case Format.FALSE:
                        callbackBool(obj, false);
                        goto Lpush;
                    default:
                        throw new UnpackException("Unknown type");
                    }

                    cur++;
                    goto Lstart;
                }
            } else {
                // data lack for deserialization
                if (used_ - cur < trail)
                    goto Labort;

                const base = cur; cur += trail - 1;  // fix current position

                final switch (state) {
                case State.FLOAT:
                    _f temp;

                    temp.i = load32To!uint(buffer_[base..base + trail]);
                    callbackFloat(obj, temp.f);
                    goto Lpush;
                case State.DOUBLE:
                    _d temp;

                    temp.i = load64To!ulong(buffer_[base..base + trail]);
                    callbackFloat(obj, temp.f);
                    goto Lpush;
                case State.REAL:
                    const expb = base + ulong.sizeof;

                    version (NonX86)
                    {
                        CustomFloat!80 temp;

                        const frac = load64To!ulong (buffer_[base..expb]);
                        const exp  = load16To!ushort(buffer_[expb..expb + ushort.sizeof]);

                        temp.significand = frac;
                        temp.exponent    = exp & 0x7fff;
                        temp.sign        = exp & 0x8000 ? true : false;

                        // NOTE: temp.get!real is inf on non-x86 when deserialized value is larger than double.max.
                        callbackFloat(obj, temp.get!real);
                    }
                    else
                    {
                        _r temp;

                        temp.fraction = load64To!(typeof(temp.fraction))(buffer_[base..expb]);
                        temp.exponent = load16To!(typeof(temp.exponent))(buffer_[expb..expb + temp.exponent.sizeof]);

                        callbackFloat(obj, temp.f);
                    }

                    goto Lpush;
                case State.UINT8:
                    callbackUInt(obj, buffer_[base]);
                    goto Lpush;
                case State.UINT16:
                    callbackUInt(obj, load16To!ushort(buffer_[base..base + trail]));
                    goto Lpush;
                case State.UINT32:
                    callbackUInt(obj, load32To!uint(buffer_[base..base + trail]));
                    goto Lpush;
                case State.UINT64:
                    callbackUInt(obj, load64To!ulong(buffer_[base..base + trail]));
                    goto Lpush;
                case State.INT8:
                    callbackInt(obj, cast(byte)buffer_[base]);
                    goto Lpush;
                case State.INT16:
                    callbackInt(obj, load16To!short(buffer_[base..base + trail]));
                    goto Lpush;
                case State.INT32:
                    callbackInt(obj, load32To!int(buffer_[base..base + trail]));
                    goto Lpush;
                case State.INT64:
                    callbackInt(obj, load64To!long(buffer_[base..base + trail]));
                    goto Lpush;
                case State.RAW: Lraw:
                    hasRaw_ = true;
                    callbackRaw(obj, buffer_[base..base + trail]);
                    goto Lpush;

                case State.EXT_DATA: Lext:
                    hasRaw_ = true;
                    obj.via.ext.type = buffer_[base];
                    callbackExt(obj, buffer_[base+1..base+trail]);
                    goto Lpush;
                case State.EXT8:
                    trail = buffer_[base] + 1;
                    if (trail == 0)
                        goto Lext;
                    state = State.EXT_DATA;
                    cur++;
                    goto Lstart;
                case State.EXT16:
                    trail = load16To!size_t(buffer_[base..base+trail]) + 1;
                    if (trail == 0)
                        goto Lext;
                    state = State.EXT_DATA;
                    cur++;
                    goto Lstart;
                case State.EXT32:
                    trail = load32To!size_t(buffer_[base..base+trail]) + 1;
                    if (trail == 0)
                        goto Lext;
                    state = State.EXT_DATA;
                    cur++;
                    goto Lstart;

                case State.STR8, State.BIN8:
                    trail = buffer_[base];
                    if (trail == 0)
                        goto Lraw;
                    state = State.RAW;
                    cur++;
                    goto Lstart;
                case State.RAW16, State.BIN16:
                    trail = load16To!size_t(buffer_[base..base + trail]);
                    if (trail == 0)
                        goto Lraw;
                    state = State.RAW;
                    cur++;
                    goto Lstart;
                case State.RAW32, State.BIN32:
                    trail = load32To!size_t(buffer_[base..base + trail]);
                    if (trail == 0)
                        goto Lraw;
                    state = State.RAW;
                    cur++;
                    goto Lstart;
                case State.ARRAY16:
                    size_t length = load16To!size_t(buffer_[base..base + trail]);
                    if (length == 0) {
                        callbackArray(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Array"(ContainerElement.ARRAY_ITEM, length);
                        state = State.HEADER;
                        cur++;
                        continue;
                    }
                case State.ARRAY36:
                    size_t length = load32To!size_t(buffer_[base..base + trail]);
                    if (length == 0) {
                        callbackArray(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Array"(ContainerElement.ARRAY_ITEM, length);
                        state = State.HEADER;
                        cur++;
                        continue;
                    }
                case State.MAP16:
                    size_t length = load16To!size_t(buffer_[base..base + trail]);
                    if (length == 0) {
                        callbackMap(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Map"(ContainerElement.MAP_KEY, length);
                        state = State.HEADER;
                        cur++;
                        continue;
                    }
                case State.MAP32:
                    size_t length = load32To!size_t(buffer_[base..base + trail]);
                    if (length == 0) {
                        callbackMap(obj, 0);
                        goto Lpush;
                    } else {
                        startContainer!"Map"(ContainerElement.MAP_KEY, length);
                        state = State.HEADER;
                        cur++;
                        continue;
                    }
                case State.HEADER:
                    break;
                }
            }

          Lpush:
            if (top == 0)
                goto Lfinish;

            auto container = &(*stack)[top - 1];

            final switch (container.type) {
            case ContainerElement.ARRAY_ITEM:
                container.value.via.array ~= obj;
                if (--container.count == 0) {
                    obj = container.value;
                    top--;
                    goto Lpush;
                }
                break;
            case ContainerElement.MAP_KEY:
                container.key  = obj;
                container.type = ContainerElement.MAP_VALUE;
                break;
            case ContainerElement.MAP_VALUE:
                container.value.via.map[container.key] = obj;
                if (--container.count == 0) {
                    obj = container.value;
                    top--;
                    goto Lpush;
                }
                container.type = ContainerElement.MAP_KEY;
            }

            state = State.HEADER;
            cur++;
        } while (cur < used_);

        goto Labort;

      Lfinish:
        (*stack)[0].value = obj;
        ret = true;
        cur++;
        goto Lend;

      Labort:
        ret = false;

      Lend:
        context_.state = state;
        context_.trail = trail;
        context_.top   = top;
        parsed_       += cur - offset_;
        offset_        = cur;

        return ret;
    }


    /**
     * supports foreach. One loop provides $(D Unpacked) object contains execute() result.
     * This is convenient in case that $(D MessagePack) values are continuous.
     */
    int opApply(scope int delegate(ref Unpacked) dg)
    {
        int result;

        while (execute()) {
            auto unpackedResult = Unpacked(context_.stack[0].value);
            result = dg(unpackedResult);
            if (result)
                break;

            clear();
        }

        return result;
    }


  private:
    /*
     * initializes internal stack environment.
     */
    @safe
    nothrow void initializeContext()
    {
        context_.state        = State.HEADER;
        context_.trail        = 0;
        context_.top          = 0;
        context_.stack.length = 1;
    }
}


unittest
{
    import msgpack.packer;

    {
        // serialize
        mixin DefinePacker;

        packer.packArray(null, true, 1, -2, "Hi!", [1], [1:1], double.max, ExtValue(7, [1,2,3,4]));

        // deserialize
        auto unpacker = StreamingUnpacker(packer.stream.data); unpacker.execute();
        auto unpacked = unpacker.purge();

        // Range test
        foreach (unused; 0..2) {
            uint i;

            foreach (obj; unpacked)
                i++;

            assert(i == unpacked.via.array.length);
        }

        auto result = unpacked.via.array;

        assert(result[0].type          == Value.Type.nil);
        assert(result[1].via.boolean   == true);
        assert(result[2].via.uinteger  == 1);
        assert(result[3].via.integer   == -2);
        assert(result[4].via.raw       == [72, 105, 33]);
        assert(result[5].as!(int[])    == [1]);
        assert(result[6].as!(int[int]) == [1:1]);
        assert(result[7].as!(double)   == double.max);
        assert(result[8].as!(ExtValue) == ExtValue(7, [1,2,3,4]));
    }

    // Test many combinations of EXT
    {
        mixin DefinePacker;

        alias Lengths = TypeTuple!(0, 1, 2, 3, 4, 5, 8, 15, 16, 31,
                                   255, 256, 2^^16, 2^^32);

        // Initialize a bunch of ExtValues and pack them
        ExtValue[Lengths.length] values;
        foreach (I, L; Lengths)
            values[I] = ExtValue(7, new ubyte[](L));
        packer.pack(values);

        auto unpacker = StreamingUnpacker(packer.stream.data); unpacker.execute();
        auto unpacked = unpacker.purge();

        // Compare unpacked values to originals
        size_t i = 0;
        foreach (deserialized; unpacked)
            assert(deserialized == values[i++]);
    }
}


private:
@trusted:


/**
 * Sets value type and value.
 *
 * Params:
 *  value = the value to set
 *  number = the content to set
 */
void callbackUInt(ref Value value, ulong number)
{
    value.type         = Value.Type.unsigned;
    value.via.uinteger = number;
}


/// ditto
void callbackInt(ref Value value, long number)
{
    value.type        = Value.Type.signed;
    value.via.integer = number;
}


/// ditto
void callbackFloat(ref Value value, real number)
{
    value.type         = Value.Type.floating;
    value.via.floating = number;
}


/// ditto
void callbackRaw(ref Value value, ubyte[] raw)
{
    value.type    = Value.Type.raw;
    value.via.raw = raw;
}

/// ditto
void callbackExt(ref Value value, ubyte[] raw)
{
    value.type    = Value.Type.ext;
    value.via.ext.data = raw;
}

/// ditto
void callbackArray(ref Value value, size_t length)
{
    value.type = Value.Type.array;
    value.via.array.length = 0;
    value.via.array.reserve(length);
}


/// ditto
void callbackMap(ref Value value, lazy size_t length)
{
    value.type    = Value.Type.map;
    value.via.map = null;  // clears previous result avoiding 'Access Violation'
}


/// ditto
void callbackNil(ref Value value)
{
    value.type = Value.Type.nil;
}


/// ditto
void callbackBool(ref Value value, bool boolean)
{
    value.type        = Value.Type.boolean;
    value.via.boolean = boolean;
}


unittest
{
    Value value;

    // Unsigned integer
    callbackUInt(value, uint.max);
    assert(value.type         == Value.Type.unsigned);
    assert(value.via.uinteger == uint.max);

    // Signed integer
    callbackInt(value, int.min);
    assert(value.type        == Value.Type.signed);
    assert(value.via.integer == int.min);

    // Floating point
    callbackFloat(value, real.max);
    assert(value.type         == Value.Type.floating);
    assert(value.via.floating == real.max);

    // Raw
    callbackRaw(value, cast(ubyte[])[1]);
    assert(value.type    == Value.Type.raw);
    assert(value.via.raw == cast(ubyte[])[1]);

    // Array
    Value[] array; array.reserve(16);

    callbackArray(value, 16);
    assert(value.type               == Value.Type.array);
    assert(value.via.array.capacity == array.capacity);

    // Map
    Value[Value] map;

    callbackMap(value, 16);
    assert(value.type    == Value.Type.map);
    assert(value.via.map == null);

    // NIL
    callbackNil(value);
    assert(value.type == Value.Type.nil);

    // Bool
    callbackBool(value, true);
    assert(value.type        == Value.Type.boolean);
    assert(value.via.boolean == true);
}
