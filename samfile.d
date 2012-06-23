module samfile;

import std.stdio;
import std.array;
import std.string;

import alignment;
import samheader;
import reference;
import sam.sam_alignment;

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

    auto alignments() @property {
        struct Result {
            this(ref File file, ref SamHeader header) {
                _header = header;
                _line_range = file.byLine();
            }
            
            bool empty() @property {
                return _line_range.empty;
            }
            
            void popFront() @property {
                _line_range.popFront();
            }

            Alignment front() @property {
                return parseAlignmentLine(cast(string)_line_range.front, _header);
            }

            alias File.ByLine!(char, char) LineRange;
            LineRange _line_range;
            SamHeader _header;

            char[] buffer;
        }
        return Result(_file, _header);
    }
private:

    File _file;

    ulong _header_end_offset;

    SamHeader _header = void;
    ReferenceSequenceInfo[] _reference_sequences;

    void _initializeStream(string filename) {
        _file = File(filename); 

        char[] _buffer;

        auto header = appender!(char[])(); 
        while (!_file.eof()) {
            _header_end_offset = _file.tell();
            auto read = _file.readln(_buffer);
            auto line = _buffer[0 .. read];

            if (line.length > 0 && line[0] == '@') {
                header.put(line);
            } else {
                if (line.length > 0) {
                    _file.seek(_header_end_offset);
                    break;
                }
            }
        }

        _header = SamHeader(cast(string)(header.data));

        foreach (sq; _header.sq_lines) {
            ReferenceSequenceInfo seq;
            seq.name = sq.sequence_name;
            seq.length = sq.sequence_length;
            _reference_sequences ~= seq;
        }
    }
}
