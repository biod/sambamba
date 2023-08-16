// Written in the D programming language.

/**
 * Compares direct conversion with stream.
 */

import std.datetime;
import std.stdio;
import std.typecons;

import msgpack;


void main()
{
    // tuple
    auto test1 = tuple(new int[](100), "MessagePack!", [1:2.0, 3:4.0, 5:6.0, 7:8.0]);
    auto data1 = pack(test1);

    // stream
    void s1()
    {
        auto result = unpack(data1).as!(typeof(test1));
    }

    // direct conversion
    void d1()
    {
        typeof(test1) result;
        unpack(data1, result);
    }

    // array
    auto test2 = new int[](1000);
    auto data2 = pack(test2);

    // stream
    void s2()
    {
        auto result = unpack(data2).as!(typeof(test2));
    }

    // direct conversion
    void d2()
    {
        typeof(test2) result;
        unpack(data2, result);
    }

    auto times = benchmark!(s1, d1, s2, d2)(1000);
    writeln("Stream(Tuple):", times[0].msecs);
    writeln("Direct(Tuple):", times[1].msecs);
    writeln("Stream(Array):", times[2].msecs);
    writeln("Direct(Array):", times[3].msecs);
}
