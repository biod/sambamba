/**
  D standard library formatting functions turned out to be
  too slow for big data processing, while standard C functions
  are very fast. This module contains set of functions for
  building strings in memory and outputting them into a file.

  For outputting to file, FILE* pointers are supported.
  For building strings in memory, provide char[] array or char*
  pointer when you're sure that amount of preallocated memory
  is enough to store string representation.

  In case char* is used, it is passed by reference, and the pointer
  is updated during string building.

  Use pointer version when it allows you to get better performance,
  but remember that it's quite dangerous.
 */
module utils.format;

import std.c.stdio;
import std.c.stdlib;
import std.string;
import std.traits;
import std.array;

/// Used for building a string in a buffer or for outputting to stream.
size_t append(Args...)(FILE* stream, string format, Args args)
{
    auto _format = toStringz(format);
    return fprintf(stream, _format, args);
}

/// ditto
size_t append(Args...)(ref char* stream, string format, Args args) {
    auto _format = toStringz(format);
    auto sz = sprintf(stream, _format, args);
    stream += sz;
    return sz;
}

/// ditto
size_t append(Args...)(ref char[] stream, string format, Args args) 
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
    return count;
}

/// Put char into a stream
size_t putcharacter(FILE* stream, char c)
{
    fputc(c, stream);    
    return 1;
}

/// Append char to a buffer
size_t putcharacter(ref char[] stream, char c)
{
    stream ~= c;
    return 1;
}

/// ditto
size_t putcharacter(ref char* stream, char c) {
    *stream++ = c;
    return 1;
}

/// Append string to output stream.
size_t putstring(T)(FILE* stream, T[] s) 
    if (is(Unqual!T == char))
{
    fwrite(s.ptr, s.length, char.sizeof, stream);
    return s.length;
}

/// Append string to a buffer
size_t putstring(T)(ref char[] stream, T[] s) 
    if (is(Unqual!T == char)) 
{
    stream ~= s;
    return s.length;
}

/// ditto
size_t putstring(T)(ref char* stream, T[] s) 
    if (is(Unqual!T == char)) 
{
    stream[0 .. s.length] = s;
    stream += s.length;
    return s.length;
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

        strreverse(str,wstr-1);

        return wstr - str;
    }
}


/// Put integer number into a stream
size_t putinteger(T)(FILE* stream, T number)
    if (isIntegral!T)
{
    char[64] buf;
    auto len = itoa(number, buf.ptr);
    fwrite(buf.ptr, len, char.sizeof, stream);
    return len;
}

/// Add string representation of an integer to a buffer
size_t putinteger(T)(ref char[] stream, T number)
    if (isIntegral!T)
{
    char[64] buf;
    auto len = itoa(number, buf.ptr);
    stream ~= buf[0 .. len];
    return len;
}

/// ditto
size_t putinteger(T)(ref char* stream, T number) {
    auto len = itoa(number, stream);
    stream += len;
    return len;
}

unittest {
    char[] buf;
    append(buf, "%d%g", 1, 2.4);
    assert(cast(string)(buf) == "12.4");

    append(buf, "%s", toStringz("k"));
    assert(cast(string)(buf) == "12.4k");

    auto str = "m";
    append(buf, "%.*s", str.length, str.ptr);
    assert(cast(string)(buf) == "12.4km");

    append(buf, "%c%c", '/', 'h');
    assert(cast(string)(buf) == "12.4km/h");

    ushort k = 5;
    append(buf, "%d", k);
    assert(cast(char)(buf[$-1]) == '5');

    buf.length = 0;
    putstring(buf, "tes");
    putcharacter(buf, 't');

    uint z = 345;
    append(buf, "%d", z);

    assert(cast(string)(buf) == "test345");

    buf.length = 0;
    putinteger(buf, 25);
    assert(cast(string)(buf) == "25");
    putinteger(buf, -31);
    assert(cast(string)(buf) == "25-31");

    char* s = cast(char*)malloc(100);
    scope(exit) free(s);

    char* p = s;
    putstring(p, "123");
    putinteger(p, 456);
    putcharacter(p, '7');
    append(p, "%g", 8.9);
    assert(s[0 .. p - s] == "12345678.9");
}
