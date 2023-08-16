// Written in the D programming language.

/**
 * Converting to/from JSON usage
 */

import std.stdio;

import msgpack;
import std.json;


void main()
{
    struct Simple
    {
        int a = 5;
        string b = "hello";
        double c = 3.14;
        char d = '!';
        @nonPacked string e = "world";
        bool f = true;
    }
    auto simple = Simple();
    Value val = simple.pack().unpack();
    writeln(val.toJSONValue());
    val = simple.pack!true().unpack();
    writeln(val.toJSONValue());

    string jsonString = `[30, 30.5, true, "hello", "40"]`;
    val = parseJSON(jsonString).fromJSONValue();
    assert(val.pack() !is null);
    assert(val.type == Value.Type.array);
    foreach (v; val.via.array)
    {
        writeln(v.type);
    }
}
