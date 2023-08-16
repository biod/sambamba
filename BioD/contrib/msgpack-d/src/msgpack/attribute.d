module msgpack.attribute;

import std.typetuple; // will use std.meta
import std.traits;


/**
 * Attribute for specifying non pack/unpack field.
 * This is an alternative approach of MessagePackable mixin.
 *
 * Example:
 * -----
 * struct S
 * {
 *     int num;
 *     // Packer/Unpacker ignores this field;
 *     @nonPacked string str;
 * }
 * -----
 */
struct nonPacked {}


package template isPackedField(alias field)
{
    enum isPackedField = (staticIndexOf!(nonPacked, __traits(getAttributes, field)) == -1) && (!isSomeFunction!(typeof(field)));
}


/**
 * Attribute for specifying serialize/deserialize proxy for pack/unpack field.
 * This is an alternative approach of registerPackHandler/registerUnpackHandler.
 *
 * Example:
 * -----
 * struct Proxy
 * {
 *     import std.conv;
 *     static void serialize(ref Packer p, ref int val) { p.pack(to!string(val)); }
 *     static void deserialize(ref Unpacker u, ref int val) { string tmp; u.unpack(tmp); val = to!int(tmp); }
 * }
 * struct S
 * {
 *     // The Packer/Unpacker proxy handler is applied this field.
 *     @serializedAs!Proxy int num;
 *     string str;
 * }
 * -----
 */
struct serializedAs(T){}

package enum bool isSerializedAs(A) = is(A : serializedAs!T, T);
package alias getSerializedAs(T : serializedAs!Proxy, Proxy) = Proxy;
package alias ProxyList(alias value) = staticMap!(getSerializedAs, Filter!(isSerializedAs, __traits(getAttributes, value)));
package template isSerializedAs(alias value)
{
    static if ( __traits(compiles, __traits(getAttributes, value)) ) {
        enum bool isSerializedAs = ProxyList!value.length > 0;
    } else {
        enum bool isSerializedAs = false;
    }
}
package template getSerializedAs(alias value)
{
    private alias _list = ProxyList!value;
    static assert(_list.length <= 1, `Only single serialization proxy is allowed`);
    alias getSerializedAs = _list[0];
}
package template hasSerializedAs(alias value)
{
    private enum _listLength = ProxyList!value.length;
    static assert(_listLength <= 1, `Only single serialization proxy is allowed`);
    enum bool hasSerializedAs = _listLength == 1;
}

unittest
{
    import msgpack.packer, msgpack.unpacker;
    struct Proxy
    {
        static void serialize(ref Packer p, ref int value) {}
        static void deserialize(ref Unpacker u, ref int value) {}
    }
    struct A
    {
        @serializedAs!Proxy int a;
        @(42) int b;
        @(42) @serializedAs!Proxy int c;
    }
    static assert(is(getSerializedAs!(A.a) == Proxy));
    static assert(isSerializedAs!(__traits(getAttributes, A.a)[0]));
    static assert(hasSerializedAs!(A.a));
    static assert(!hasSerializedAs!(A.b));
    static assert(hasSerializedAs!(A.c));
}
