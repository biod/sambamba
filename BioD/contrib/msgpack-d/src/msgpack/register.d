module msgpack.register;

import msgpack.packer;
import msgpack.unpacker;

import std.array;


/**
 * Register a serialization handler for $(D_PARAM T) type
 *
 * Example:
 * -----
 * registerPackHandler!(Foo, fooPackHandler);
 * -----
 */
void registerPackHandler(T, alias Handler, Stream = Appender!(ubyte[]))()
{
    PackerImpl!(Stream).registerHandler!(T, Handler);
}


/**
 * Register a deserialization handler for $(D_PARAM T) type
 *
 * Example:
 * -----
 * registerUnackHandler!(Foo, fooUnackHandler);
 * -----
 */
void registerUnpackHandler(T, alias Handler)()
{
    Unpacker.registerHandler!(T, Handler);
}


/**
 * Register derived class for (de)serialization
 *
 * Example:
 * -----
 * registerClass!(DerivedClass);
 * -----
 */
void registerClass(T, Stream = Appender!(ubyte[]))()
{
    PackerImpl!(Stream).register!(T);
    Unpacker.register!(T);
}
