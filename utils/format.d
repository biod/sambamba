/**
  D standard library formatting functions turned out to be
  too slow for big data processing, while standard C functions
  are very fast. This module contains set of functions for
  building strings in memory and outputting them into a file.
 */
module utils.format;

import std.c.stdio;
import std.c.stdlib;
import std.string;

/// Used for building a string in a buffer or for outputting to stream.
void append(Args...)(FILE* stream, string format, Args args)
{
    auto _format = toStringz(format);
    fprintf(stream, _format, args);
}

/// ditto
void append(Args...)(ref char[] stream, string format, Args args) 
{
    char[1024] buffer;
    int count;

    auto f = toStringz(format);
    auto p = buffer.ptr;
    auto psize = buffer.length;
    for (;;)
    {
        version(Win32)
        {
            count = _snprintf(p,psize,f,args);
            if (count != -1)
                break;
            psize *= 2;
            p = cast(char *) alloca(psize);
        }
        version(Posix)
        {
            count = snprintf(p,psize,f,args);
            if (count == -1)
                psize *= 2;
            else if (count >= psize)
                psize = count + 1;
            else
                break;
            p = cast(char *) alloca(psize);
        }
    }

    if (stream.capacity < stream.length + count) {
        stream.reserve(stream.capacity * 2);
    }
    stream ~= p[0 .. count].dup;
}

/// Put a char into a stream
void putcharacter(FILE* stream, char c)
{
    fputc(c, stream);    
}

/// ditto
void putcharacter(ref char[] stream, char c)
{
    stream ~= c;
}

/// Appends string to output stream or buffer.
void putstring(Out)(ref Out stream, string s) {
    append(stream, "%.*s", s.length, s.ptr);
}

unittest {
    char[] buf;
    append(buf, "%d%g", 1, 2.4);
    assert(buf == "12.4");

    append(buf, "%s", toStringz("k"));
    assert(buf == "12.4k");

    auto str = "m";
    append(buf, "%.*s", str.length, str.ptr);
    assert(buf == "12.4km");

    append(buf, "%c%c", '/', 'h');
    assert(buf == "12.4km/h");

    ushort k = 5;
    append(buf, "%d", k);
    assert(buf[$-1] == '5');

    buf = "".dup;
    putstring(buf, "tes");
    putchar(buf, 't');

    uint z = 345;
    append(buf, "%d", z);

    assert(buf == "test345");
}
