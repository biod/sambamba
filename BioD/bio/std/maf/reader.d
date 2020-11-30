/*
    This file is part of BioD.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>
    Copyright (C) 2018    Pjotr Prins <pjotr.prins@thebird.nl>

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
module bio.std.maf.reader;

import bio.std.maf.block;
import bio.std.maf.parser;

import std.array;
import std.string;
import std.stdio;
import std.algorithm;
import std.range.interfaces;
import std.string;

///
struct MafBlockRange {
    private {

        File _f;
        bool _empty;
        MafBlock _front;

        void skipHeader() {
          auto _lines = _f.byLine();
          if (!_lines.empty)
            if (_lines.front.startsWith("##maf"))
                _lines.popFront();
        }
    }

    this(string fn) {
        _f = File(fn);
        // _lines = _f.byLine(KeepTerminator.yes);
        skipHeader();
        popFront();
    }

    ///
    bool empty() @property const {
        return _empty;
    }

    ///
    MafBlock front() @property {
        return _front;
    }

    ///
    void popFront() {
        auto block_data = Appender!(char[])();
        auto _lines = _f.byLine();
        while (!_lines.empty && !_lines.front.chomp().empty) {
            block_data.put(_lines.front.dup);
            _lines.popFront();
        }
        if (block_data.data.empty) {
            _empty = true;
        } else {
            _front = parseMafBlock(cast(string)(block_data.data));
            if (!_lines.empty)
                _lines.popFront();
        }
    }
}


///
class MafReader {

    private string _fn;

    ///
    this(string filename) {
        _fn = filename;
    }

    ///
    string filename() @property const {
        return _fn;
    }

    ///
    MafBlockRange blocks() @property {
        return MafBlockRange(_fn);
    }
}
