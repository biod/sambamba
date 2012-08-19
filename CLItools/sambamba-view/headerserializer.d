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
module headerserializer;

import samheader;
import std.json;
import std.conv;
import std.array;
import std.c.stdio;
import utils.format;
import utils.msgpack;

private {
    interface IHeaderSerializer {
        void writeln(SamHeader header);
    }

    final class HeaderSamSerializer : IHeaderSerializer {
        void writeln(SamHeader header) {
            serialize(header, stdout);
        }
    }

    final class HeaderJsonSerializer : IHeaderSerializer {

        void writeln(SamHeader header) {
            
            static JSONValue jv(T)(T value) {
                JSONValue v;
                static if(is(T == string)) {
                    v.str = value;
                    v.type = JSON_TYPE.STRING;
                } else static if(is(T == uint)) {
                    v.integer = value;
                    v.type = JSON_TYPE.INTEGER;
                }
                return v;
            }

            JSONValue result;
            result.type = JSON_TYPE.OBJECT;

            result.object["format_version"] = jv(header.format_version);

            if (header.sorting_order != SortingOrder.unknown) {
                result.object["sorting_order"] = jv(to!string(header.sorting_order));
            }

            JSONValue tmp;
            tmp.type = JSON_TYPE.ARRAY;
            tmp.array = new JSONValue[header.sequences.length];
            size_t i = 0;
            foreach (line; header.sequences) {
                JSONValue sq;
                sq.type = JSON_TYPE.OBJECT;
                sq.object["sequence_name"] = jv(line.name);
                sq.object["sequence_length"] = jv(line.length);
                sq.object["assembly"] = jv(line.assembly);
                sq.object["md5"] = jv(line.md5);
                sq.object["species"] = jv(line.species);
                sq.object["uri"] = jv(line.uri);
                tmp.array[i++] = sq;
            }
            result.object["sq_lines"] = tmp;

            tmp.array.length = header.read_groups.length;
            i = 0;
            foreach (line; header.read_groups) {
                JSONValue sq;
                sq.type = JSON_TYPE.OBJECT;
                sq.object["identifier"] = jv(line.identifier);
                sq.object["sequencing_center"] = jv(line.sequencing_center);
                sq.object["description"] = jv(line.description);
                sq.object["date"] = jv(line.date);
                sq.object["flow_order"] = jv(line.flow_order);
                sq.object["key_sequence"] = jv(line.key_sequence);
                sq.object["library"] = jv(line.library);
                sq.object["programs"] = jv(line.programs);
                sq.object["predicted_insert_size"] = jv(line.predicted_insert_size);
                sq.object["platform"] = jv(line.platform);
                sq.object["platform_unit"] = jv(line.platform_unit);
                sq.object["sample"] = jv(line.sample);
                tmp.array[i++] = sq;
            }
            result.object["rg_lines"] = tmp;

            tmp.array.length = header.programs.length;
            i = 0;
            foreach (line; header.programs) {
                JSONValue sq;
                sq.type = JSON_TYPE.OBJECT;
                sq.object["identifier"] = jv(line.identifier);
                sq.object["program_name"] = jv(line.name);
                sq.object["command_line"] = jv(line.command_line);
                sq.object["previous_program"] = jv(line.previous_program);
                sq.object["program_version"] = jv(line.program_version);
                tmp.array[i++] = sq;
            }
            result.object["pg_lines"] = tmp;

            putstring(stdout, toJSON(&result));
            putcharacter(stdout, '\n');
        }
    }
}

final class HeaderMsgpackSerializer : IHeaderSerializer {
    void writeln(SamHeader header) {
        auto packer = packer(Appender!(ubyte[])());

        void packArrayOf(T)(T objs) {
            packer.beginArray(objs.length);
            foreach (obj; objs)
                packer.pack(obj);
        }

        packer.pack(header.format_version);
        packer.pack(to!string(header.sorting_order));
        packArrayOf(header.sequences);
        packArrayOf(header.read_groups);
        packArrayOf(header.programs);

        fwrite(packer.stream.data.ptr, packer.stream.data.length, ubyte.sizeof, stdout);
    }
}

final class HeaderSerializer {
    private IHeaderSerializer _serializer;

    this(string format) {
        switch(format) {
            case "sam":
            case "bam":
                _serializer = new HeaderSamSerializer(); break;
            case "json":
                _serializer = new HeaderJsonSerializer(); break;
            case "msgpack":
                _serializer = new HeaderMsgpackSerializer(); break;
            default:
                throw new Exception("unknown format for serialization: '" ~ format ~ 
                                    "' (expected 'sam', 'bam', 'json', or 'msgpack')");
        }
    }

    void writeln(SamHeader header) {
        _serializer.writeln(header);
    }
}
