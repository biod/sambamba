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
module bio.std.hts.bam.referenceinfo;

import contrib.undead.stream;
import std.exception;
import std.array;

/**
  Stores basic information about reference sequence.
 */
struct ReferenceSequenceInfo {
    private {
        string _name;
        int _length;
    }

    /// Reference sequence name
    /// (null byte is guaranteed to follow the returned slice)
    string name() @property const nothrow {
        return _name[0 .. $ - 1];
    }
   
    /// Reference sequence length
    int length() @property const {
        return _length;
    }

    ///
    this(string name, int length) {
        _name = name ~ '\0';
        _length = length;
    }

    /// Constructs the structure from input stream
    this(ref Stream stream) {
        int l_name; // length of the reference name plus one
        stream.read(l_name);
        _name = cast(string)stream.readString(l_name); // keep '\0' at the end
        stream.read(_length);
    }
}
