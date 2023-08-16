// Written in the D programming language.

/**
 * Attribute usage
 */

import msgpack;

struct Hoge
{
    string f1;
    @nonPacked int f2;
}

void main()
{
    Hoge hoge = Hoge("hoge", 10);
    Hoge fuga;

    unpack(pack(hoge), fuga);
    assert(hoge.f1 == fuga.f1);
    assert(hoge.f2 != fuga.f2);
    assert(fuga.f2 == int.init);
}
