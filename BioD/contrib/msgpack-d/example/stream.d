// Written in the D programming language.

/**
 * Stream deserialization usage
 */

import std.array;
import std.concurrency;
import std.exception;
import std.stdio;

import msgpack;


void deserializer()
{
    auto unpacker = StreamingUnpacker(cast(ubyte[])null);
    bool endLoop;

    while (true) {
        receive((immutable(ubyte)[] data) { unpacker.feed(data); },
                (bool end) { endLoop = end; });

        if (endLoop)
            break;

        while (unpacker.execute()) {
            auto unpacked = unpacker.purge();
            writeln("Type:  ", unpacked.type);
            writeln("Value: ", unpacked.as!(string));
        }

        if (unpacker.size >= 100)
            throw new Exception("Too large!");
    }
}


void main()
{
    string message = "Hell";
    foreach (i; 0..93)  // Throws Exception if 94
        message ~= 'o';

    auto packed = pack(message);
    auto data   = packed.assumeUnique();
    auto tid    = spawn(&deserializer);

    while (!data.empty) {
        auto limit = data.length >= 10 ? 10 : data.length;

        tid.send(data[0..limit]);
        data = data[limit..$];
    }

    tid.send(true);
}
