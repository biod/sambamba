/*
    This file is part of BioD.
    Copyright (C) 2012-2015    Artem Tarasov <lomereiter@gmail.com>

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
module bio.std.hts.bam.writer;

import bio.std.hts.bam.referenceinfo;
import bio.std.hts.sam.header;
import bio.std.hts.bam.constants;
import bio.std.hts.bam.bai.indexing;
import bio.std.hts.bam.read;
import bio.std.hts.bam.readrange;
import bio.core.bgzf.outputstream;
import bio.core.bgzf.virtualoffset;
import bio.core.utils.stream;

import std.parallelism;
import std.exception;
import contrib.undead.stream;
import std.traits;
import std.system;
import std.algorithm;
import std.array;
import std.bitmanip;

/** Class for outputting BAM.
    $(BR)
    Compresses BGZF blocks in parallel.
    Tries to write reads so that they don't cross BGZF block borders.
    $(BR)
    Usage is very simple, see example below.

    Example:
    --------------------------------------
    import bio.std.hts.bam.writer, bio.std.hts.bam.reader;
    ...
    auto src = new BamReader("in.bam");
    auto dst = new BamWriter("out.bam", 9); // maximal compression
    scope (exit) dst.finish();              // close the stream at exit
    dst.writeSamHeader(src.header);         // copy header and reference sequence info
    dst.writeReferenceSequenceInfo(src.reference_sequences);
    foreach (read; src.reads) {
        if (read.mapping_quality > 10)      // skip low-quality reads
            dst.writeRecord(read);
    }
    --------------------------------------
    */
final class BamWriter {

    /// Creates new BAM writer outputting to file or $(I stream).
    /// Automatically writes BAM magic number (4 bytes).
    ///
    /// Params:
    ///     compression_level  = compression level, must be in range -1..9
    ///     task_pool          = task pool to use for parallel compression
    ///     buffer_size        = size of BgzfOutputStream buffer
    this(contrib.undead.stream.Stream stream,
         int compression_level=-1,
         std.parallelism.TaskPool task_pool=std.parallelism.taskPool,
         size_t buffer_size=0)
    {
        _stream = new BgzfOutputStream(stream, compression_level,
                                       task_pool, buffer_size);
        _stream.setWriteHandler((ubyte[] uncompressed, ubyte[] compressed) {
            _bytes_written += compressed.length;
        });

        writeString(BAM_MAGIC);
    }

    /// ditto
    this(string filename,
         int compression_level=-1,
         std.parallelism.TaskPool task_pool=std.parallelism.taskPool)
    {
        _filename = filename;
        auto filestream = new bio.core.utils.stream.File(filename, "wb+");
        this(filestream, compression_level, task_pool);
    }

    /// Can be called right after the stream constructor, only once
    void setFilename(string output_filename) {
        enforce(_filename is null, "Can't set output filename twice");
        _filename = output_filename;
    }

    /// By default, the writer attempts to automatically create index
    /// when writing coordinate-sorted files. If this behaviour is not
    /// desired, it can be switched off before writing SAM header.
    void disableAutoIndexCreation() {
        _disable_index_creation = true;
    }

    package void writeByteArray(const(ubyte[]) array) {
        _stream.writeExact(array.ptr, array.length);
    }

    package void writeString(string str) {
        writeByteArray(cast(ubyte[])str);
    }

    package void writeInteger(T)(T integer) if (isIntegral!T)
    {
        ubyte[T.sizeof] buf = nativeToLittleEndian(integer);
        _stream.writeExact(buf.ptr, buf.length);
    }

    private {
        size_t _bytes_written;
        bool _create_index = false;
        bool _disable_index_creation = false;
        bool _record_writing_mode = false;
        string _filename;

        IndexBuilder _index_builder;

        VirtualOffset _start_vo, _end_vo;
        ubyte[] _pending_read_data_buf;
        ubyte[] _pending_read_data;

        void appendReadData(ubyte[] data) {
            auto required = _pending_read_data.length + data.length;
            if (_pending_read_data_buf.length < required)
                _pending_read_data_buf.length = max(required, _pending_read_data_buf.length * 2);
            _pending_read_data_buf[_pending_read_data.length .. $][0 .. data.length] = data;
            _pending_read_data = _pending_read_data_buf[0 .. required];
        }

        size_t _len;

        void indexBlock(ubyte[] uncompressed, ubyte[] compressed) {
            ushort inner_offset = 0;

            void indexRead(ubyte[] data) {
                auto read = BamRead(data);
                if (uncompressed.length > 0) {
                    _end_vo = VirtualOffset(_bytes_written, inner_offset);
                } else {
                    _end_vo = VirtualOffset(_bytes_written + compressed.length, 0);
                }
                auto read_block = BamReadBlock(_start_vo, _end_vo, read);
                _index_builder.put(read_block);
            }

            if (_pending_read_data !is null) {
                if (uncompressed.length < _len) {
                    appendReadData(uncompressed);
                    _len -= uncompressed.length;
                    uncompressed = null;
                } else {
                    appendReadData(uncompressed[0 .. _len]);
                    uncompressed = uncompressed[_len .. $];
                    inner_offset = cast(ushort)_len;
                    indexRead(_pending_read_data);
                    _pending_read_data = null;
                }
            }

            while (uncompressed.length > 0) {
                _len = *cast(int*)(uncompressed.ptr); // assume LE...
                _start_vo = VirtualOffset(_bytes_written, inner_offset);
                if (_len + int.sizeof <= uncompressed.length) {
                    _pending_read_data = null;
                    auto read_data = uncompressed[int.sizeof .. int.sizeof + _len];
                    uncompressed = uncompressed[int.sizeof + _len .. $];
                    inner_offset += _len + int.sizeof;
                    indexRead(read_data);
                } else { // read spans multiple BGZF blocks
                    appendReadData(uncompressed[int.sizeof .. $]);
                    _len -= _pending_read_data.length;
                    break;
                }
            }
            _bytes_written += compressed.length;
        }
    }

    /// Writes SAM header. Should be called after construction.
    void writeSamHeader(bio.std.hts.sam.header.SamHeader header) {
        writeSamHeader(header.text);
    }

    /// ditto
    void writeSamHeader(string header_text) {
        _create_index = !_disable_index_creation &&
            !header_text.find("SO:coordinate").empty &&
            _filename.length >= 4 &&
            _filename[$ - 4 .. $] == ".bam";
        writeInteger(cast(int)header_text.length);
        writeString(header_text);
    }

    /// Writes reference sequence information. Should be called after
    /// dumping SAM header. Writer will store this array to use later for
    /// resolving read reference IDs to names.
    ///
    /// Flushes current BGZF block.
    void writeReferenceSequenceInfo(const(bio.std.hts.bam.referenceinfo.ReferenceSequenceInfo)[] reference_sequences)
    {
        _reference_sequences = reference_sequences;

        auto n_refs = cast(int)reference_sequences.length;
        writeInteger(n_refs);
        foreach (sequence; reference_sequences) {
            writeInteger(cast(int)(sequence.name.length + 1));
            writeString(sequence.name);
            writeInteger(cast(ubyte)'\0');
            writeInteger(cast(int)sequence.length);
        }

        if (_create_index) {
            auto index = new bio.core.utils.stream.File(_filename ~ ".bai", "wb+");
            _index_builder = IndexBuilder(index, n_refs);
            _index_builder.check_bins = true;
        }

        _stream.flushCurrentBlock();
    }

    private void indexingWriteHandler(ubyte[] uncomp, ubyte[] comp) {
        indexBlock(uncomp, comp);
    }

    /// Writes BAM read. Throws exception if read reference ID is out of range.
    void writeRecord(R)(R read) {
        enforce(read.ref_id == -1 || read.ref_id < _reference_sequences.length,
                "Read reference ID is out of range");

        if (!_record_writing_mode) {
            if (_create_index) {
                _record_writing_mode = true;
                _stream.setWriteHandler(&indexingWriteHandler);
            } else {
                _stream.setWriteHandler(null);
            }
        }

        read._recalculate_bin();

        auto read_size = read.size_in_bytes;
        if (read_size + _current_size > BGZF_BLOCK_SIZE) {
            _stream.flushCurrentBlock();
            read.write(this);
            _current_size = read_size;
        } else {
            read.write(this);
            _current_size += read_size;
        }
    }

    /// Flushes current BGZF block.
    void flush() {
        _stream.flush();
    }

    /// Flushes buffer and closes output stream. Adds BAM EOF block automatically.
    void finish() {
        _stream.close();
        if (_create_index)
            _index_builder.finish();
    }

    private {
        BgzfOutputStream _stream;
        const(ReferenceSequenceInfo)[] _reference_sequences;
        size_t _current_size; // number of bytes written to the current BGZF block
    }
}
