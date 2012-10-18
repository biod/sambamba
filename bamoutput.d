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
module bamoutput;

import samheader;
import reference;
import alignment;
import alignmentrange;

import bgzfcompress;
import constants;

import std.stream;
import std.system;
import std.range;
import std.traits;
import std.algorithm;
import std.parallelism;
import std.exception;

debug {
    import std.stdio;
}

/// Writes BAM magic string at the beginning of the file.
void writeMagic(EndianStream stream) {
    stream.writeString(BAM_MAGIC);
}

/// Writes header length (int32_t) and header text to a stream.
void writeSamHeader(EndianStream stream, string header) 
in
{
    assert(stream.endian == Endian.littleEndian);
}
body
{
    stream.write(cast(int)(header.length));
    stream.writeString(header);
}

/// Writes reference sequence information to a stream.
void writeReferenceSequences(EndianStream stream, ReferenceSequenceInfo[] info) 
in
{
    assert(stream.endian == Endian.littleEndian);
}
body
{
    stream.write(cast(int)info.length);
    foreach (refseq; info) {
        stream.write(cast(int)(refseq.name.length + 1));
        stream.writeString(refseq.name);
        stream.write(cast(ubyte)'\0');
        stream.write(cast(int)(refseq.length));
    }
}

/// Writes alignment to a stream.
void writeAlignment(EndianStream stream, Alignment alignment) {
    alignment.write(stream);
}

/// Writes BAM file to a stream.
///
/// Params:
///
///     stream     =  the stream to write to
///
///     header     =  SAM header as raw text
///
///     info       =  reference sequences info
///
///     alignments =  range of alignments
///
///     compression_level =  level of compression, use 0 for uncompressed BAM.
/// 
void writeBAM(R)(Stream stream, 
                 string header, 
                 ReferenceSequenceInfo[] info,
                 R alignments,
                 int compression_level=-1,
                 TaskPool task_pool=taskPool,
                 size_t parmapbufsz=32,
                 size_t parmapwusz=1) 
    if (is(Unqual!(ElementType!R) == Alignment)
            || is(Unqual!(ElementType!R) == AlignmentBlock))
{
    // First, pack header and reference sequences.
    auto header_memory_stream = new MemoryStream();
    auto header_endian_stream = new EndianStream(header_memory_stream, Endian.littleEndian);
    writeMagic(header_endian_stream);
    writeSamHeader(header_endian_stream, header);
    writeReferenceSequences(header_endian_stream, info);

    foreach (block; std.range.chunks(header_memory_stream.data, BGZF_BLOCK_SIZE)) {
        auto bgzf_block = bgzfCompress(block, compression_level);
        stream.writeExact(bgzf_block.ptr, bgzf_block.length);
    }

    // Range of blocks, each is <= MAX_BLOCK_SIZE in size,
    // except cases where single alignment takes more than
    // one block. In this particular case, the alignment occupies
    // the whole block.
    static struct BlockRange(R) {
        this(R alignments) {
            _alignments = alignments;

            _memory_stream = new MemoryStream();
            _endian_stream = new EndianStream(_memory_stream);

            if (!_alignments.empty) {
                _current_alignment = _alignments.front;
            }

            loadNextBlock();
        }

        R _alignments;
        bool _empty = false;

        MemoryStream _memory_stream;
        EndianStream _endian_stream;

        Alignment _current_alignment;

        bool empty() @property {
            return _empty;
        }

        ubyte[] front() @property {
            return _memory_stream.data[0 .. cast(size_t)_endian_stream.position];
        }

        void popFront() {
            _endian_stream.seekSet(0);
            loadNextBlock();
        }

        void loadNextBlock() {
            if (_alignments.empty) {
                _empty = true;
                return;
            }
            while (true) {
                auto pos = _endian_stream.position;
                if (pos == 0 || (_current_alignment.size_in_bytes + pos <= BGZF_BLOCK_SIZE)) 
                {
                    writeAlignment(_endian_stream, _current_alignment);
                    _alignments.popFront();
                    if (!_alignments.empty) {
                        _current_alignment = _alignments.front;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
        }
    }

    static auto blockRange(R)(R alignments) {
        return BlockRange!R(alignments);
    }

    auto blocks = blockRange(alignments);
    auto chunked_blocks = chunkedBlockRange(blocks, BGZF_BLOCK_SIZE);

    // ugly workaround 
    // (issue 5710, cannot use delegates as parameters to non-global template)
    static ubyte[] makeBgzfCompressor(int n)(ubyte[] buf) {
        return bgzfCompress(buf, n);
    }
       
    // helper function
    static void writeAlignmentBlocks(int n, R)(R chunked_blocks, ref Stream stream,
                                               TaskPool task_pool,
                                               size_t parmapbufsz,
                                               size_t parmapwusz) {

    version(serial) {
        auto bgzf_blocks = map!(makeBgzfCompressor!n)(chunked_blocks);
    } else {
        auto bgzf_blocks = task_pool.map!(makeBgzfCompressor!n)(chunked_blocks,
                                                                parmapbufsz,
                                                                parmapwusz);
    }

        foreach (bgzf_block; bgzf_blocks) {
            stream.writeExact(bgzf_block.ptr, bgzf_block.length);
        }
    }

    switch (compression_level) {
        case -1: writeAlignmentBlocks!(-1)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 0: writeAlignmentBlocks!(0)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 1: writeAlignmentBlocks!(1)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 2: writeAlignmentBlocks!(2)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 3: writeAlignmentBlocks!(3)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 4: writeAlignmentBlocks!(4)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 5: writeAlignmentBlocks!(5)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 6: writeAlignmentBlocks!(6)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 7: writeAlignmentBlocks!(7)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 8: writeAlignmentBlocks!(8)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        case 9: writeAlignmentBlocks!(9)(chunked_blocks, stream, task_pool, parmapbufsz, parmapwusz); break;
        default: throw new Exception("compression level must be a number from -1 to 9");
    }

    // write EOF block
    stream.writeExact(BAM_EOF.ptr, BAM_EOF.length);
}

private {
    // Takes range of blocks as input and returns
    // another range where too big blocks, with length
    // more than $(D max_size) were cut into several smaller blocks
    struct ChunkedBlockRange(R) {
        R _blocks;
        ElementType!R _current_block;
        size_t _max_size;
        bool _empty = false;
        this(R blocks, size_t max_size) {
            _blocks = blocks;
            _max_size = max_size;
            _loadNextBlock();
        }

        bool empty() @property {
            return _empty;
        }

        ubyte[] front() @property {
            if (_current_block.length <= _max_size) {
                return _current_block;
            } else {
                version(serial) {
                    return _current_block[0 .. _max_size];
                } else {
                    return _current_block[0 .. _max_size].dup;
                }
            }
        }

        void popFront() {
            if (_current_block.length <= _max_size) {
                _loadNextBlock();
            } else {
                _current_block = _current_block[_max_size .. $];
            }
        }

        void _loadNextBlock() {
            if (!_blocks.empty) {
                // dup is required because map.front returns a reference
                _current_block = _blocks.front.dup;
                _blocks.popFront();
            } else {
                _empty = true;
            }
        }
    }

    auto chunkedBlockRange(R)(R blocks, size_t max_size) {
        return ChunkedBlockRange!R(blocks, max_size);
    }

}

unittest {
    import std.stdio;
    writeln("Testing chunkedBlockRange...");
    ubyte[][] blocks = [[1, 2, 3], [4, 5, 6, 7], [8, 9, 10]];
    assert(equal(chunkedBlockRange(blocks, 2), 
                 [[1, 2], [3], [4, 5], [6, 7], [8, 9], [10]]));
    assert(equal(chunkedBlockRange(blocks, 3),
                 [[1, 2, 3], [4, 5, 6], [7], [8, 9, 10]]));
    assert(equal(chunkedBlockRange(blocks, 4), blocks));
    assert(equal(chunkedBlockRange(blocks, 1), map!"[a]"(iota(1, 11))));
}
