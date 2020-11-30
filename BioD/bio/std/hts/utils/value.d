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
module bio.std.hts.utils.value;

import bio.std.hts.bam.tagvalue;
import std.exception;
import std.array;

private {
    static ubyte* arrayPointer(ref Value v) {
        // if value contains array, as a tagged union, it has the following layout:
        //
        // Value = size_t
        //         void*
        //         other stuff...
        return cast(ubyte*)(*(cast(size_t*)(&v) + 1));
    }
}

/// Emplace value at address $(D p).
/// Assumes that enough memory is allocated at that address.
/// (You can find needed amount of memory with $(D sizeInBytes) function)
void emplaceValue(ubyte* p, ref Value v) {
    enforce(!v.is_nothing, "null value can't be stored in BAM");

    auto tag = v.tag;
    auto size = tag >> 5; // element size

    if ((tag & 1) == 0) { // primitive type
        *p++ = cast(ubyte)v.bam_typeid;

        p[0 .. size] = (cast(ubyte*)(&v))[0 .. size];
    } else {
        
        auto bytes = *cast(size_t*)(&v) * (tag >> 5);

        if (v.is_string) {
            *p++ = cast(ubyte)v.bam_typeid;
            
            p[0 .. bytes] = arrayPointer(v)[0 .. bytes];
            p[bytes] = 0; // trailing zero

        } else {
            *p++ = cast(ubyte)'B';
            *p++ = cast(ubyte)v.bam_typeid;

            *(cast(uint*)p) = cast(uint)(bytes / size); // # of elements
            p += uint.sizeof;

            p[0 .. bytes] = arrayPointer(v)[0 .. bytes];
        }
    }
}

/// Put a value to an Appender!(ubyte[]) struct
void emplaceValue(ref Appender!(ubyte[]) appender, ref Value v) {
    enforce(!v.is_nothing, "null value can't be stored in BAM");

    auto tag = v.tag;
    auto size = tag >> 5; // element size

    if ((tag & 1) == 0) { // primitive type
        appender.put(cast(ubyte)v.bam_typeid);
        appender.put((cast(ubyte*)(&v))[0 .. size]);
    } else {
        
        auto bytes = *cast(size_t*)(&v) * (tag >> 5);

        if (v.is_string) {
            appender.put(cast(ubyte)v.bam_typeid);
            appender.put(arrayPointer(v)[0..bytes]); 
            appender.put(cast(ubyte)0); // trailing zero
        } else {
            appender.put(cast(ubyte)'B');
            appender.put(cast(ubyte)v.bam_typeid);
            uint number_of_elems = cast(uint)(bytes / size);
            appender.put((cast(ubyte*)(&number_of_elems))[0..4]);
            appender.put(arrayPointer(v)[0 .. bytes]);
        }
    }
}

/// Calculate size in bytes which value will consume in BAM file.
size_t sizeInBytes(ref Value v) {
    enforce(!v.is_nothing, "null value can't be stored in BAM");

    auto tag = v.tag;

    if ((tag & 1) == 0) {
        return char.sizeof + (tag >> 5); // primitive type
    } 

    auto bytes = *cast(size_t*)(&v) * (tag >> 5);

    if (v.is_string) {
        return char.sizeof + bytes + char.sizeof; // trailing zero
    } else {
        return 2 * char.sizeof + uint.sizeof + bytes;
    }
}
