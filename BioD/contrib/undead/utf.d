/**
 * Contains the obsolete functions from Phobos' `std.utf`.
 */
module contrib.undead.utf;

import std.utf;
import std.typecons;

//deprecated("Removed October 2017. Please use std.utf.encode instead.")
char[] toUTF8(return out char[4] buf, dchar c) nothrow @nogc @safe pure
{
    const sz = encode!(Yes.useReplacementDchar)(buf, c);
    return buf[0 .. sz];
}

//deprecated("Removed October 2017. Please use std.utf.encode instead.")
wchar[] toUTF16(return ref wchar[2] buf, dchar c) nothrow @nogc @safe pure
{
    const sz = encode!(Yes.useReplacementDchar)(buf, c);
    return buf[0 .. sz];
}
