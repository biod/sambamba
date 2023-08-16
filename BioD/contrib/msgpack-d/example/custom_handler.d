import msgpack;
import std.array;
import std.stdio;
import std.variant;

class A { }
class C : A 
{
    int num;
    this(int n) { num = n; }
}

void cPackHandler(ref Packer p, ref C c)
{
    writeln("Pack C: ", c.num);
    p.pack(c.num);
}

void cUnpackHandler(ref Unpacker u, ref C c)
{
    writeln("Unpack C: ", c.num);
    u.unpack(c.num);
}

void vPackHandler(ref Packer p, ref Variant v)
{
    writeln("pack Variant: ", v);
    p.pack(v.get!bool);
}

void vUnpackHandler(ref Unpacker u, ref Variant v)
{
    writeln("unpack Variant: ", v);
    bool b;
    u.unpack(b);
    v = b;
}

void main()
{
    registerPackHandler!(C, cPackHandler);
    registerUnpackHandler!(C, cUnpackHandler);
    registerPackHandler!(Variant, vPackHandler);
    registerUnpackHandler!(Variant, vUnpackHandler);

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

        Variant v = true;
        p.pack(v);

        Variant v2 = 10;
        unpack(p.stream.data, v2);
        assert(v2 == true);
    }
}
