// Written in the D programming language.

/**
 * Stream Deserializer with foreach.
 */

import std.stdio;

import msgpack;


void main()
{
    // create 3 MessagePack objects([1, 0.1L], true, "foobarbaz")
    auto test1 = pack(1, 0.1L) ~ pack(true);
    auto test2 = pack("foobarbaz");

    // split data to deserialize test
    test1 ~= test2[0..2];
    test2  = test2[2..$];

    auto unpacker = StreamingUnpacker(test1);

    foreach (unpacked; unpacker) {
        if (unpacked.type == Value.Type.array) {
            foreach (obj; unpacked) {
                switch (obj.type) {
                case Value.Type.unsigned: writeln(obj.as!(uint)); break;
                case Value.Type.floating: writeln(obj.as!(real)); break;
                default:
                    throw new Exception("Unknown type");
                }
            }
        } else {
            writeln(unpacked.as!(bool));
        }
    }

    unpacker.feed(test2);

    foreach (unpacked; unpacker)
        writeln(unpacked.as!(string));
}
