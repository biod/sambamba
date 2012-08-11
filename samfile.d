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

extern(C) size_t lseek(int fd, size_t offset, int whence);

import alignment;
import samheader;
import reference;
import sam.recordparser;

struct SamFile {

    this(string filename, TaskPool task_pool=taskPool) {
        _file = File(filename);
        _filename = filename;
        _seekable = lseek(_file.fileno, 0, 0) != ~0;
        _initializeStream(filename);
        _task_pool = task_pool;
    }

    SamHeader header() @property {
        return _header;
    }

    ReferenceSequenceInfo[] reference_sequences() @property {
        return _reference_sequences;
    }

    private alias File.ByLine!(char, char) LineRange;

    /// Alignments in SAM file. Not thread-safe.
    /// If file is seekable, starts from the beginning,
    /// otherwise range starts with the alignment first
    /// in the stream.
    auto alignments() @property {

        LineRange _lines;
        if (_seekable) {
            File _file;
            if (_filename is null) {
                this._file.seek(0);
                _file = this._file;
            } else {
                _file = File(_filename);
            }
            _lines = _file.byLine();
            auto dummy = _lines.front;
            for (long i = 0; i < _lines_to_skip && !_lines.empty; i++) {
                _lines.popFront();
            }
        } else {
            _lines = this._lines;
        }

        auto lines = map!"a.idup"(_lines);

        version (serial) {

            auto build_storage = new AlignmentBuildStorage();
            Alignment parse(string s) {
                return parseAlignmentLine(s, _header, build_storage);
            }

            return map!parse(lines);

        } else {

            static __gshared SamHeader header;
            header = _header;
            return _task_pool.map!((string s) { return parseAlignmentLine(s, header); })(lines, 1024);

        }
    }
private:

    File _file;
    string _filename;
    bool _seekable;
    long _lines_to_skip;
    LineRange _lines;
    TaskPool _task_pool;

    SamHeader _header;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream(string filename) {
        auto header = appender!(char[])(); 

        _lines = _file.byLine();

        while (!_lines.empty) {
            auto line = _lines.front;
            if (line.length > 0 && line[0] == '@') {
                header.put(line);
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
