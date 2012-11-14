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

import alignment;
import samheader;
import reference;
import sam.recordparser;

private {
    extern(C) size_t lseek(int, size_t, int);
}

struct SamFile {

    this(string filename) {
        _file = File(filename);
        _filename = filename;
        _seekable = lseek(_file.fileno, 0, 0) != ~0;
        _initializeStream();
    }

    SamHeader header() @property {
        return _header;
    }

    ReferenceSequenceInfo[] reference_sequences() @property {
        return _reference_sequences;
    }

    private alias File.ByLine!(char, char) LineRange;

    static struct SamRecordRange {
        this(LineRange lines, ref SamHeader header) {
            _header = header;
            _line_range = lines;

            _build_storage = new AlignmentBuildStorage();
            _parseNextLine();
        }
        
        bool empty() @property {
            return _empty;
        }
        
        void popFront() @property {
            _line_range.popFront();
            _parseNextLine();
        }

        Alignment front() @property {
            return _current_alignment;
        }

        private {
            void _parseNextLine() {
                if (_line_range.empty) {
                    _empty = true;
                } else {
                    _current_alignment = parseAlignmentLine(cast(string)_line_range.front.dup,
                                                            _header,
                                                            _build_storage);
                }
            }

            LineRange _line_range;
            Alignment _current_alignment;
            bool _empty;
            SamHeader _header;
            AlignmentBuildStorage _build_storage;
        }
    }

    /// Alignments in SAM file. Can be iterated only once.
    auto alignments() @property {
        
        LineRange lines = _lines;
        if (_seekable) {
            if (_filename !is null) {
                File file = File(_filename);
                lines = file.byLine();
            } else {
                _file.seek(0);
                lines = _file.byLine();
            }
            auto dummy = lines.front;
            for (int i = 0; i < _lines_to_skip; i++)
                lines.popFront();
        }

        return SamRecordRange(lines, _header);
    }
private:

    File _file;
    bool _seekable;
    string _filename;
    LineRange _lines;
    ulong _lines_to_skip;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream() {
        auto header = appender!(char[])(); 

        _lines = _file.byLine();

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
                header.put('\n');
                _lines_to_skip += 1;
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
