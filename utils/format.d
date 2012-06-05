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
import std.traits;

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

    stream ~= p[0 .. count];
}

/// Put char into a stream
void putcharacter(FILE* stream, char c)
{
    fputc(c, stream);    
}

/// Append char to a buffer
void putcharacter(ref char[] stream, char c)
{
    stream ~= c;
}

/// Append string to output stream.
void putstring(T)(FILE* stream, T[] s) 
    if (is(Unqual!T == char))
{
    fwrite(s.ptr, s.length, char.sizeof, stream);
}

/// Append string to a buffer
void putstring(T)(ref char[] stream, T[] s) 
    if (is(Unqual!T == char)) 
{
    stream ~= s;
}

private {
    /// Reverses closed interval [begin .. end]
    void strreverse(char* begin, char* end)
    {
        char aux;
        while (end > begin)
            aux = *end, *end-- = *begin, *begin++ = aux;
    }

    /// Prints $(D value) at the address where $(D str) points to.
    /// Returns number of characters written.
    auto itoa(T)(T value, char* str)
    {
        char* wstr=str;

        static if (isSigned!T) { 
            uint uvalue = (value < 0) ? -value : value;
        } else {
            uint uvalue = value;
        }

        do {
            *wstr++ = cast(char)(48 + (uvalue % 10)); 
        } while (uvalue /= 10);

        static if (isSigned!T) {
            if (value < 0) *wstr++ = '-';
        }

        // Reverse string
        strreverse(str,wstr-1);

        return wstr - str;
    }
}


/// Put integer number into a stream
void putinteger(T)(FILE* stream, T number)
    if (isIntegral!T)
{
    char[64] buf;
    auto len = itoa(number, buf.ptr);
    fwrite(buf.ptr, len, char.sizeof, stream);
}

/// Add string representation of an integer to a buffer
void putinteger(T)(ref char[] stream, T number)
    if (isIntegral!T)
{
    char[64] buf;
    auto len = itoa(number, buf.ptr);
    stream ~= buf[0..len];
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
    putcharacter(buf, 't');

    uint z = 345;
    append(buf, "%d", z);

    assert(buf == "test345");

    buf = [];
    putinteger(buf, 25);
    assert(buf == "25");
    putinteger(buf, -31);
    assert(buf == "25-31");
}
