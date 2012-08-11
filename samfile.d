/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module samfile;

import std.stdio;
import std.array;
import std.string;
import std.algorithm;
import std.parallelism;

import alignment;
import samheader;
import reference;
import sam.recordparser;

struct SamFile {

    this(string filename) {
        _initializeStream(filename);
    }

    SamHeader header() @property {
        return _header;
    }

    ReferenceSequenceInfo[] reference_sequences() @property {
        return _reference_sequences;
    }

    private alias File.ByLine!(char, char) LineRange;

    /// Alignments in SAM file. Can be iterated only once.
    auto alignments() @property {
        auto build_storage = new AlignmentBuildStorage();
        auto lines = map!"a.idup"(_lines);

        Alignment parse(string s) {
            return parseAlignmentLine(s, _header, build_storage);
        }

        return map!parse(lines); // TODO: parallelize
    }
private:

    File _file;
    LineRange _lines;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream(string filename) {
        _file = File(filename); 

        char[] _buffer;

        auto header = appender!(char[])(); 

        _lines = _file.byLine();

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
                _lines.popFront();
            } else {
                break;
            }
        }

        _header = new SamHeader(cast(string)(header.data));

        _reference_sequences = new ReferenceSequenceInfo[_header.sequences.length];
        foreach (sq; _header.sequences) {
            ReferenceSequenceInfo seq;
            seq.name = sq.name;
            seq.length = sq.length;
            _reference_sequences[_header.getSequenceIndex(seq.name)] = seq;
        }
    }
}
