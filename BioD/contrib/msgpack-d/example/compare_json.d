// Written in the D programming language.

/**
 * Compares std.json
 */

import std.datetime;
import std.json;
import std.stdio;

import msgpack;


void main()
{
    JSONValue jsonObj = parseJSON(`[12, "foo", true, 0.23, {"1":1}, [1, 2]]`);

    void f1()
    {
        parseJSON(toJSON(jsonObj));
    }

    Value mpObj = unpack(pack(12, "foo", true, 0.23, ["1":1], [1, 2]));

    void f2()
    {
        unpack(pack(mpObj));
    }

    auto times = benchmark!(f1, f2)(10000);
    writeln("JSON:    ", times[0].msecs);
    writeln("Msgpack: ", times[1].msecs);
}
