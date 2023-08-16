import std.conv;
import std.stdio;
import std.datetime;
import msgpack;

struct A {
    int x;
}

struct Foo
{
    A[string] a;
}

void main()
{
    Foo foo;
    foreach (a; 'a' .. 'z')
        foreach (b; 'a' .. 'z')
            foreach (c; 'a' .. 'z')
                    foo.a[to!string(a) ~ to!string(b) ~ to!string(c)] = A();

    auto sw = StopWatch(AutoStart.yes);
    ubyte[] data = msgpack.pack(foo);
    writeln(sw.peek.usecs);
    
    auto sw2 = StopWatch(AutoStart.yes);
    Foo foo1 = msgpack.unpack(data).as!(typeof(foo));
    writeln(sw2.peek.usecs);

    assert(foo == foo1);

    Foo foo2;
    auto sw3 = StopWatch(AutoStart.yes);
    msgpack.unpack(data, foo2);
    writeln(sw3.peek.usecs);

    assert(foo == foo2);
}
