/*
   This file is part of BioD.
   Copyright (C) 2016   George Githinji <biorelated@gmail.com>
   Copyright (C) 2018   Emilio Palumbo <emiliopalumbo@gmail.com> 
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

module bio.std.file.fai;

import std.stdio;
import std.range;
import std.algorithm;
import std.string;
import std.conv;
import std.file;
import std.path;

struct FaiRecord {
    string header, lineTerm;
    ulong seqLen, lineLen, offset;
    @property ulong lineOffset() {
        return lineLen + lineTerm.length;
    }

    string toString() {
        return format("%s\t%s\t%s\t%s\t%s", header, seqLen, offset, lineLen, lineOffset);
    }
    unittest {
        auto rec = FaiRecord("chr2", "\n", 10, 50, 4);
        assert(rec.toString() == "chr2\t10\t4\t50\t51");
        rec.lineTerm = "\r\n";
        assert(rec.toString() == "chr2\t10\t4\t50\t52");
    }

    this(string str) {
        auto res = str.split("\t");
        header ~= res[0];
        seqLen = to!ulong(res[1]);
        offset = to!ulong(res[2]);
        lineLen = to!ulong(res[3]);
        lineTerm = (to!ulong(res[4])-lineLen) == 1 ? "\n" : "\r\n";
    }
    unittest {
        auto s = "chr2\t10\t4\t50\t51";
        assert(FaiRecord(s).toString() == s);
    }

    this(string header, string lineTerm, ulong seqLen, ulong lineLen, ulong offset) {
        this.header = header;
        this.seqLen = seqLen;
        this.offset = offset;
        this.lineLen = lineLen;
        this.lineTerm = lineTerm;
    }
    unittest {
        assert(FaiRecord("chr2", "\n", 10, 50, 4).toString() == "chr2\t10\t4\t50\t51");
    }
}

auto readFai(string filename) {
    File f = File(filename, "r");
    return f.byLineCopy()
            .map!(x => FaiRecord(x));
}
unittest {
    auto faiString = "chr2\t10\t4\t50\t51";
    auto testIndex = tempDir.buildPath("test.fa.fai");
    scope(exit) testIndex.remove;
    File(testIndex, "w").writeln(faiString);
    auto recs = readFai(testIndex).array;
    assert(recs.length == 1);
    assert(is(typeof(recs[0])==FaiRecord));
    assert(recs[0].toString() == faiString);
}

auto makeIndex(T)(T records) {   
    FaiRecord[string] index;
    foreach (record; records) {
        index[record.header] = record;
    }
    index.rehash;
    return index;
}
unittest {
    auto records = to!(FaiRecord[])(["chr2\t10\t4\t50\t51"]);
    auto i = makeIndex(records);
    assert( i.length == 1);
    assert( "chr2" in i);
    assert( i["chr2"] ==  FaiRecord("chr2\t10\t4\t50\t51"));
}

auto buildFai(string filename) {

    File f = File(filename, "r");
    FaiRecord[] records; 
    string lineTerm = f.byLine(KeepTerminator.yes).take(1).front.endsWith("\r\n") ? "\r\n" : "\n";
    f.seek(0);
    ulong offset;
    foreach(line; f.byLine(KeepTerminator.no, lineTerm)) {
        offset+= line.length + lineTerm.length;
        if ( line.startsWith(">") ) {
            records~=FaiRecord();
            records[$-1].lineTerm = lineTerm;
            records[$-1].header ~= line.split(" ").front[1..$];
            records[$-1].offset = offset;
        } else {
            if ( records[$-1].lineLen == 0 ) {
                records[$-1].lineLen = line.length;
            }
            records[$-1].seqLen += line.length;
        }
    }

    return records;
}

unittest {
    auto testFa = tempDir.buildPath("test.fa");
    scope(exit) testFa.remove;
    File(testFa, "w").writeln(q"(
        >chr1
        acgtgagtgc
        >chr2
        acgtgagtgcacgtgagtgcacgtgagtgc
        acgtgagtgcacgtgagtgc
    )".outdent().strip());
    auto recs = buildFai(testFa).array;
    assert(recs.length == 2);
    assert(recs.all!(x => is(typeof(x)==FaiRecord)));
    assert(recs[0].toString() == "chr1\t10\t6\t10\t11");
    assert(recs[1].toString() == "chr2\t50\t23\t30\t31");
}