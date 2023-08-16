[![CI](https://github.com/msgpack/msgpack-d/actions/workflows/d.yml/badge.svg)](https://github.com/msgpack/msgpack-d/actions/workflows/d.yml)

# MessagePack for D

MessagePack is a binary-based JSON-like serialization library.

MessagePack for D is a pure D implementation of MessagePack.

# Features

* Small size and High performance
* Zero copy serialization / deserialization
* Streaming deserializer for non-contiguous IO situation
* Supports D features (Ranges, Tuples, real type)

Note: The `real` type is only supported in D.
Don't use the `real` type when communicating with other programming languages.
Note that `Unpacker` will raise an exception if a loss of precision occurs.

## Current Limitations

* No circular references support
* If you want to use the LDC compiler, you need at least version 0.15.2 beta2

# Install

Use dub to add it as a dependency:

```sh
% dub install msgpack-d
```

# Usage

Example code can be found in the `example` directory.

The documentation can be found [here](http://msgpack.github.io/msgpack-d/)

## pack / unpack

msgpack-d is very simple to use. Use `pack` for serialization, and `unpack` for deserialization:

```D
import std.file;
import msgpack;

struct S { int x; float y; string z; }

void main()
{
    S input = S(10, 25.5, "message");

    // serialize data
    ubyte[] inData = pack(input);

    // write data to a file
    write("file.dat", inData);

    // read data from a file
    ubyte[] outData = cast(ubyte[])read("file.dat");

    // unserialize the data
    S target = outData.unpack!S();

    // verify data is the same
    assert(target.x == input.x);
    assert(target.y == input.y);
    assert(target.z == input.z);
}
```

### Feature: Skip serialization/deserialization of a specific field.

Use the `@nonPacked` attribute:

```d
struct User
{
    string name;
    @nonPacked int level;  // pack / unpack will ignore the 'level' field
}
```

### Feature: Use your own serialization/deserialization routines for custom class and struct types.

msgpack-d provides the functions `registerPackHandler` / `registerUnpackHandler` to allow you
to use custom routines during the serialization or deserialization of user-defined class and struct types.
This feature is especially useful when serializing a derived class object when that object is statically
typed as a base class object.

For example:

```d
class Document { }
class XmlDocument : Document
{
    this() { }
    this(string name) { this.name = name; }
    string name;
}

void xmlPackHandler(ref Packer p, ref XmlDocument xml)
{
    p.pack(xml.name);
}

void xmlUnpackHandler(ref Unpacker u, ref XmlDocument xml)
{
    u.unpack(xml.name);
}

void main()
{
    /// Register the 'xmlPackHandler' and 'xmlUnpackHandler' routines for
    /// XmlDocument object instances.
    registerPackHandler!(XmlDocument, xmlPackHandler);
    registerUnpackHandler!(XmlDocument, xmlUnpackHandler);

    /// Now we can serialize/deserialize XmlDocument object instances via a
    /// base class reference.
    Document doc = new XmlDocument("test.xml");
    auto data = pack(doc);
    XmlDocument xml = unpack!XmlDocument(data);
    assert(xml.name == "test.xml");  // xml.name is "test.xml"
}
```

In addition, here is also a method using `@serializedAs` attribute:

```d
import std.datetime: Clock, SysTime;
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

void main()
{
    /// Now we can serialize/deserialize LogData
    LogData[] logs;
    logs ~= LogData("MessagePack is nice!");
    auto data = pack(logs);
    LogData[] datas = unpack!(LogData[])(data);
    assert(datas[0].timestamp.toString() == datas[0].timestamp.toString());
}
```

## The PackerImpl / Unpacker / StreamingUnpacker types

These types are used by the `pack` and `unpack` functions.

See the documentation of [PackerImpl](http://msgpack.github.io/msgpack-d/#PackerImpl), [Unpacker](http://msgpack.github.io/msgpack-d/#Unpacker) and [StreamingUnpacker](http://msgpack.github.io/msgpack-d/#StreamingUnpacker) for more details.

# Links

* [The MessagePack Project](http://msgpack.org/)

  The official MessagePack protocol website.

* [msgpack-d's issue tracker](https://github.com/msgpack/msgpack-d/issues)

  Use this issue tracker to review and file bugs in msgpack-d.

* [MessagePack's Github](http://github.com/msgpack/)

  Other language bindings and implementations of the msgpack protocol can be found here.

# Copyright

    Copyright (c) 2010- Masahiro Nakagawa

# License

Distributed under the [Boost Software License, Version 1.0](http://www.boost.org/users/license.html).
