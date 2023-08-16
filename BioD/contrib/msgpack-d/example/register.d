import msgpack;
import std.stdio;

class A
{
    string str = "foo";
}

class C : A 
{
    int num;
    this(int n) { num = n; }
}

void main()
{
    registerClass!(C);

    {
        Packer p;
        A c = new C(1000);
        p.pack(c);

        A c2 = new C(5);
        unpack(p.stream.data, c2);
        assert(1000 == (cast(C)c2).num);
    }
    {
        Packer p;
        C c = new C(1000);
        p.pack(c);

        C c2 = new C(5);
        unpack(p.stream.data, c2);
        assert(1000 == (cast(C)c2).num);
    }
}
