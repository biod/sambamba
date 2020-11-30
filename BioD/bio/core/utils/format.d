/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
/**
  $(P This module provides fast formatting functions.)

  $(P Each function has two overloads:
    $(UL
        $(LI $(D ref char*) - in this case, function starts
            writing at the location, and updates the pointer.
            No checks are done, it's user's responsibility that this is safe.)
        $(LI $(D scope void delegate(const(char)[])) - formatted data
            is passed to the delegate for further processing.)))
 */
module bio.core.utils.format;

import core.stdc.stdio;
import core.stdc.stdlib;
static import core.stdc.string;
import std.string;
import std.traits;
import std.array;
import std.math;

///
template isSomeSink(T) {
    static if (__traits(compiles, T.init("string")))//T == void delegate(const(char)[])))
        enum isSomeSink = true;
    else static if (is(T == char*))
        enum isSomeSink = true;
    else
        enum isSomeSink = false;
}

private {
    // Reverses closed interval [begin .. end]
    void strreverse(char* begin, char* end)
    {
        char aux;
        while (end > begin)
            aux = *end, *end-- = *begin, *begin++ = aux;
    }

    // Prints $(D value) at the address where $(D str) points to.
    // Returns number of characters written.
    size_t itoa(T)(T value, char* str)
    {
        char* wstr=str;

        static if (isSigned!T) {
            ulong uvalue = (value < 0) ? -cast(int)(value) : value;
        } else {
            ulong uvalue = value;
        }

        do {
            *wstr++ = cast(char)(48 + (uvalue % 10));
        } while (uvalue /= 10);

        static if (isSigned!T) {
            if (value < 0) *wstr++ = '-';
        }

        strreverse(str,wstr-1);

        return wstr - str;
    }
}


private {
    void writeFloat(T)(ref char* sink, T number)
        if (isFloatingPoint!T)
    {
        char[4] format;
        format[0] = '%';
        format[1] = 'g';
        format[2] = '\0';
        sink += sprintf(sink, format.ptr, number);
    }

    void writeFloat(T)(scope void delegate(const(char)[]) sink, T number)
        if (isFloatingPoint!T)
    {
        char[1024] buffer = void;
        int count;

        auto p = buffer.ptr;
        auto psize = buffer.length;
        for (;;)
        {
            version(Win32)
            {
                count = _snprintf(p,psize,"%g", cast(double)number);
                if (count != -1)
                    break;
                psize *= 2;
                p = cast(char *) alloca(psize);
            }
            version(Posix)
            {
                count = snprintf(p,psize,"%g", cast(double)number);
                if (count == -1)
                    psize *= 2;
                else if (count >= psize)
                    psize = count + 1;
                else
                    break;
                p = cast(char *) alloca(psize);
            }
        }

        sink(p[0 .. count]);
    }

    void writeInteger(T)(ref char* sink, T integer)
        if (isIntegral!T)
    {
        sink += itoa(integer, sink);
    }

    void writeInteger(T)(scope void delegate(const(char)[]) sink, T integer)
        if (isIntegral!T)
    {
        char[32] buf = void;
        auto len = itoa(integer, buf.ptr);
        sink(buf[0 .. len]);
    }

    void writeChar(T)(ref char* sink, T c)
        if (isSomeChar!T)
    {
        *sink++ = c;
    }

    void writeChar(T)(scope void delegate(const(char)[]) sink, T c)
        if (isSomeChar!T)
    {
        sink((&c)[0 .. 1]);
    }

    void writeString(T)(ref char* sink, T s)
        if (isSomeString!T)
    {
        auto str = cast(const(char)[])s;
        core.stdc.string.memcpy(sink, str.ptr, str.length);
        sink += str.length;
    }

    void writeString(T)(scope void delegate(const(char)[]) sink, T s)
        if (isSomeString!T)
    {
        sink(cast(const(char)[])s);
    }

    void writeImpl(Sink, T)(auto ref Sink sink, T value)
        if (isSomeSink!Sink)
    {
        static if (isIntegral!T)
            writeInteger(sink, value);
        else static if (isFloatingPoint!T)
            writeFloat(sink, value);
        else static if (isSomeChar!T)
            writeChar(sink, value);
        else static if (isSomeString!T)
            writeString(sink, value);
        else static assert(false,
                    "only integers, floats, chars and strings are supported");
    }

    // -------------------- JSON output utils ----------------------------------

    // JSON doesn't support NaN and +/- infinity.
    // Therefore the approach taken here is to represent
    // infinity as 1.0e+1024, and NaN as null.
    void writeFloatJson(Sink, T)(auto ref Sink sink, T value)
        if (isFloatingPoint!T)
    {
        if (isFinite(value)) {
            sink.write(value);
        } else {
            if (value == float.infinity) {
                sink.write("1.0e+1024");
            } else if (value == -float.infinity) {
                sink.write("-1.0e+1024");
            } else if (isNaN(value)) {
                sink.write("null");
            } else {
                assert(0);
            }
        }
    }

    immutable char[256] specialCharacterTable = [
    /*   0-15  */    0,0,  0,0,0,0,0,0, 'b','t','n',0, 'f','r',0,  0,
    /*  16-31  */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,   0,  0,0,  0,
    /*  32-47  */    0,0,'"',0,0,0,0,0,   0,  0,  0,0,   0,  0,0,  0,
    /*  48-63  */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,   0,  0,0,'/',
    /*  64-79  */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,   0,  0,0,  0,
    /*  80-95  */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,'\\',  0,0,  0,
    /*  96-111 */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,   0,  0,0,  0,
    /* 112-127 */    0,0,  0,0,0,0,0,0,   0,  0,  0,0,   0,  0,0,  0,

                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
    ];

    void writeStringJson(Sink, T)(auto ref Sink sink, T s)
        if (isSomeString!T)
    {
        sink.write('"');
        foreach (char c; s) {
            auto sc = specialCharacterTable[cast(ubyte)c];
            if (sc == 0) {
                sink.write(c);
            } else {
                sink.write('\\');
                sink.write(sc);
            }
        }
        sink.write('"');
    }

    void writeCharJson(Sink, T)(auto ref Sink sink, T c)
        if (isSomeChar!T)
    {
        sink.writeStringJson((&c)[0 .. 1]);
    }

    void writeArrayJson(Sink, T)(auto ref Sink sink, T array)
        if (isArray!T && __traits(compiles, sink.writeJson(array[0])))
    {
        if (array.length == 0) {
            sink.write("[]");
            return;
        }

        sink.write('[');
        foreach (elem; array[0 .. $ - 1]) {
            sink.writeJson(elem);
            sink.write(',');
        }
        sink.writeJson(array[$ - 1]);
        sink.write(']');
    }

    void writeJsonImpl(Sink, T)(auto ref Sink sink, T value)
        if (isSomeSink!Sink)
    {
        static if (isIntegral!T)
            writeInteger(sink, value);
        else static if (isFloatingPoint!T)
            writeFloatJson(sink, value);
        else static if (isSomeChar!T)
            writeCharJson(sink, value);
        else static if (isSomeString!T)
            writeStringJson(sink, value);
        else static if (isArray!T && __traits(compiles, sink.writeJsonImpl(value[0])))
            writeArrayJson(sink, value);
        else static assert(false,
                    "only numbers, chars, strings and arrays are supported");
    }
}

///
void write(T)(ref char* sink, T value) { writeImpl(sink, value); }
///
void write(T)(scope void delegate(const(char)[]) sink, T value) { writeImpl(sink, value); }

///
void writeArray(Sink, T, U)(auto ref Sink sink, T array, U delimiter)
    if (isSomeSink!Sink && isArray!T && (isSomeChar!U || isSomeString!U) &&
        __traits(compiles, sink.write(array[0])))
{
    if (array.length == 0)
        return;

    foreach (elem; array[0 .. $ - 1]) {
        sink.write(elem);
        sink.write(delimiter);
    }
    sink.write(array[$ - 1]);
}

/// Supports numbers, strings, and arrays. No dictionary - because D doesn't have a good one.
void writeJson(T)(ref char* sink, T value) { writeJsonImpl(sink, value); }
/// ditto
void writeJson(T)(scope void delegate(const(char)[]) sink, T value) { writeJsonImpl(sink, value); }
