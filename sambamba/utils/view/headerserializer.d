/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.utils.view.headerserializer;

import bio.sam.header;
import std.conv;
import std.array;
import std.stdio;
import bio.bam.thirdparty.msgpack;

private {
    abstract class IHeaderSerializer {
        private File _f;
        this(File f) { _f = f; }
        void writeln(SamHeader header);
    }

    final class HeaderSamSerializer : IHeaderSerializer {
        this(File f) { super(f); }
        override void writeln(SamHeader header) {
            _f.writef("%s", header);
        }
    }

    final class HeaderJsonSerializer : IHeaderSerializer {
        this(File f) { super(f); }
        override void writeln(SamHeader header) {
            _f.writef("%j", header);
        }
    }

    final class HeaderMsgpackSerializer : IHeaderSerializer {
        this(File f) { super(f); }
        override void writeln(SamHeader header) {
            auto packer = packer(Appender!(ubyte[])());
            packer.pack(header);
            fwrite(packer.stream.data.ptr, packer.stream.data.length, ubyte.sizeof, _f.getFP());
        }
    }
}

final class HeaderSerializer {
    private IHeaderSerializer _serializer;

    this(File f, string format) {
        switch(format) {
            case "sam":
            case "bam":
                _serializer = new HeaderSamSerializer(f); break;
            case "json":
                _serializer = new HeaderJsonSerializer(f); break;
            case "msgpack":
                _serializer = new HeaderMsgpackSerializer(f); break;
            default:
                throw new Exception("unknown format for serialization: '" ~ format ~ 
                                    "' (expected 'sam', 'bam', 'json', or 'msgpack')");
        }
    }

    void writeln(SamHeader header) {
        _serializer.writeln(header);
    }
}
