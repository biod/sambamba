module bio.std.sff.writer;

import bio.std.sff.constants;
import bio.std.sff.utils.roundup;

import bio.core.utils.stream;
import contrib.undead.stream;
import std.system;

/// Class for outputting SFF files
class SffWriter {

    /// Create new writer.
    this(string filename, string flow_order, string key_sequence) 
    {
        _filename = filename;
        _flow_order = flow_order;
        _key_seq = key_sequence;

        auto f = new bio.core.utils.stream.File(filename, "wb+");
        auto stream = new BufferedStream(f, 1024576);
        _endian_stream = new EndianStream(stream, Endian.bigEndian);

        writeHeader();
    }

    /// Flow order
    string flow_order() @property const {
        return _flow_order;
    }

    /// Key sequence
    string key_sequence() @property const {
        return _key_seq;
    }

    /// Add a read to the end of file
    void append(R)(R sff_read) {
        // compute read_header_length
        ushort exact_read_header_length = cast(ushort)(16 + sff_read.name.length);
        ushort read_header_length = roundup(exact_read_header_length);

        _endian_stream.write(read_header_length);
        _endian_stream.write(cast(ushort)sff_read.name.length);
        _endian_stream.write(cast(uint)sff_read.bases.length);
        _endian_stream.write(sff_read.clip_qual_left);
        _endian_stream.write(sff_read.clip_qual_right);
        _endian_stream.write(sff_read.clip_adapter_left);
        _endian_stream.write(sff_read.clip_adapter_right);
        _endian_stream.writeExact(sff_read.name.ptr, sff_read.name.length);
        for (size_t i = 0; i < read_header_length - exact_read_header_length; ++i)
            _endian_stream.write(cast(ubyte)0);

        for (size_t i = 0; i < _flow_order.length; ++i)
            _endian_stream.write(sff_read.flowgram_values[i]);

        auto n_bases = sff_read.bases.length;
        _endian_stream.writeExact(sff_read.flow_index_per_base.ptr, n_bases);
        _endian_stream.writeExact(sff_read.bases.ptr, n_bases);
        _endian_stream.writeExact(sff_read.quality_scores.ptr, n_bases);

        auto k = 2 * _flow_order.length + 3 * n_bases;
        auto padding = roundup(k) - k;
        
        for (size_t i = 0; i < padding; ++i)
            _endian_stream.write(cast(ubyte)0);

        ++_n_reads;
    }

    /// Flush all buffers and update number of reads in the file header
    void finish() {
        updateNumberOfReads();
        _endian_stream.close();
    }

    private {
        string _filename;
        string _flow_order;
        string _key_seq;
        Stream _endian_stream;

        uint _n_reads;

        ushort _exact_header_len() @property const {
            return cast(ushort)(31 + _flow_order.length + _key_seq.length);
        }

        ushort _header_len() @property const {
            return roundup(_exact_header_len);
        }

        void writeHeader() {
            _endian_stream.write(SFF_MAGIC);
            _endian_stream.writeExact(SFF_VERSION.ptr, 4);
            _endian_stream.write(0UL);
            _endian_stream.write(0U);
            _endian_stream.write(_n_reads);
            _endian_stream.write(_header_len);
            _endian_stream.write(cast(ushort)_key_seq.length);
            _endian_stream.write(cast(ushort)_flow_order.length);
            _endian_stream.write(cast(ubyte)1);
            _endian_stream.writeExact(_flow_order.ptr, _flow_order.length);
            _endian_stream.writeExact(_key_seq.ptr, _key_seq.length);
            for (size_t i = 0; i < _header_len - _exact_header_len; ++i)
                _endian_stream.write(cast(ubyte)0);
        }

        void updateNumberOfReads() {
            auto old_pos = _endian_stream.position;
            _endian_stream.position = 20;
            _endian_stream.write(_n_reads);
            _endian_stream.position = old_pos;
        }
    }
}
