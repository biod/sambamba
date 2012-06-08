module bamoutput;

import samheader;
import reference;
import alignment;

import bgzfcompress;
import constants;

import std.stream;
import std.system;
import std.range;
import std.algorithm;

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
///     compression_level =  level of compression, use 0 for uncompressed BAM;
///                          default value is 9 - maximum compression.
/// 
void writeBAM(R)(Stream stream, 
                 string header, 
                 ReferenceSequenceInfo[] info,
                 R alignments,
                 int compression_level=9)
    if (is(ElementType!R == Alignment))
{
    // First, pack header and reference sequences.
    auto header_memory_stream = new MemoryStream();
    auto header_endian_stream = new EndianStream(header_memory_stream, Endian.littleEndian);
    writeMagic(header_endian_stream);
    writeSamHeader(header_endian_stream, header);
    writeReferenceSequences(header_endian_stream, info);

    auto header_block = bgzfCompress(header_memory_stream.data, compression_level);
    stream.writeExact(header_block.ptr, header_block.length);

    // OK, now alignments

    // Range of blocks, each is <= MAX_BLOCK_SIZE in size,
    // except cases where single alignment takes more than
    // one block. In this particular case, the alignment occupies
    // the whole block.
    struct BlockRange {
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
            return _memory_stream.data[0 .. _endian_stream.position];
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

    auto blocks = BlockRange(alignments);
    foreach(block; blocks) {
        while (block.length > 0) {
            auto to_write = min(block.length, BGZF_BLOCK_SIZE);

            debug {
                writeln("compressing block of length ", to_write);
            }

            auto bgzf_block = bgzfCompress(block[0 .. to_write], compression_level);

            stream.writeExact(bgzf_block.ptr, bgzf_block.length);

            block = block[to_write .. $];
        }
    }

    // write EOF block
    stream.writeExact(BAM_EOF.ptr, BAM_EOF.length);
}
