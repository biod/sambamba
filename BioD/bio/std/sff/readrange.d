module bio.std.sff.readrange;

import bio.std.sff.read;
import bio.std.sff.index;
import bio.core.utils.switchendianness;

import std.algorithm;
import contrib.undead.stream;
import std.system;
import std.array;

private {
    
    // GC used in D is quite bad at allocating lots of objects in a tight loop.
    // The following is a simple way to reduce the number of allocations.

    ubyte[] current_chunk;
    size_t used;

    size_t chunk_size = 65_536;

    static this() {
        current_chunk = uninitializedArray!(ubyte[])(chunk_size);
        used = 0;
    }

    T[] allocateArray(T : T[])(size_t size) {
        size_t new_used = used + size * T.sizeof;
        if (new_used > chunk_size) {
            new_used = size * T.sizeof;
            if (new_used > chunk_size)
                chunk_size = new_used;

            current_chunk = uninitializedArray!(ubyte[])(chunk_size);
            used = new_used;
            return cast(T[])current_chunk[0 .. used];
        } else {
            auto old_used = used;
            used = new_used;
            return cast(T[])current_chunk[old_used .. used];
        }
    }
}

struct SffReadRange {
    this(Stream stream, 
         ushort number_of_flows_per_read,
         IndexLocation index_location)
    {
        _stream = stream;
        _n_flows = number_of_flows_per_read;
        _index_loc = index_location;

        _fetchNextRead();
    }

    private {
        Stream _stream;
        ushort _n_flows;
        IndexLocation _index_loc;

        bool _empty;
        SffRead _read;

        void _fetchNextRead() {
            if (_stream.position == _index_loc.offset)
                _stream.seekCur(_index_loc.length);

            if (_stream.eof) {
                _empty = true;
            } else {
                _read.file_offset = _stream.position;
                // determine how many bytes to read
                ushort read_header_length = void;
                ushort name_length = void;
                uint number_of_bases = void;
                
                _stream.read(read_header_length);
                _stream.read(name_length);
                _stream.read(number_of_bases);
                _stream.read(_read.clip_qual_left);
                _stream.read(_read.clip_qual_right);
                _stream.read(_read.clip_adapter_left);
                _stream.read(_read.clip_adapter_right);

                char[] name = allocateArray!(char[])(name_length);
                _stream.readExact(name.ptr, name_length);
                _stream.seekCur(read_header_length - 16 - name_length);
                _read.name = cast(string)name;

                size_t _data_length = _n_flows * ushort.sizeof + 3 * number_of_bases;

                _read.flowgram_values = allocateArray!(ushort[])(_n_flows);
                _stream.readExact(_read.flowgram_values.ptr, _n_flows * ushort.sizeof);

                if (std.system.endian != Endian.bigEndian) {
                    for (size_t i = 0; i < _n_flows; ++i) {
                        switchEndianness(_read.flowgram_values.ptr + i, ushort.sizeof);
                    }
                }

                _read.flow_index_per_base = allocateArray!(ubyte[])(number_of_bases);
                _stream.readExact(_read.flow_index_per_base.ptr, number_of_bases);

                _read.bases = allocateArray!(char[])(number_of_bases);
                _stream.readExact(_read.bases.ptr, number_of_bases);

                _read.quality_scores = allocateArray!(ubyte[])(number_of_bases);
                _stream.readExact(_read.quality_scores.ptr, number_of_bases);

                if (_data_length % 8 > 0)
                    _stream.seekCur(8 - (_data_length % 8));
            }
        }
    }

    bool empty() @property const {
        return _empty;
    }

    SffRead front() @property {
        return _read;
    }

    void popFront() {
        _fetchNextRead();
    }
}
