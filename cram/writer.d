module cram.writer;

import undead.stream;
import std.string;

import cram.htslib;
import cram.exception;

import bio.bam.read;
import bio.bam.referenceinfo;
import bio.sam.header;

extern(C) {
    int fai_build(const(char)* fn);
}

string fastaIndexFilename(string fasta_fn) {
    import std.file : exists, isFile;
    auto fn = fasta_fn ~ ".fai";
    if (fn.exists && fn.isFile)
        return fn;
    if (fai_build(toStringz(fasta_fn)) < 0)
        throw new Exception("Failed to build fasta index for " ~ fasta_fn);
    return fn;
}

final class CramWriter {
    private {
        cram_fd* _fd;
        string _fn;
        size_t _n_threads;

        string _reference_fn;
        string _fasta_idx_fn;
    }

    this(string filename, size_t n_threads) {
        _fn = filename;
        _n_threads = n_threads;

        _fd = cram_open(toStringz(filename), "wb");
        if (_fd is null)
            throw new CramException("Can't open file for writing: " ~ filename);

        if (_n_threads > 1)
            cram_set_option(_fd, cram_option.CRAM_OPT_NTHREADS, _n_threads);
    }

    void setFastaFilename(string ref_fn) {
        _reference_fn = ref_fn ~ "\0";
        _fasta_idx_fn = fastaIndexFilename(ref_fn);
    }

    void writeSamHeader(SamHeader header) {
        writeSamHeader(header.text);
    }

    void writeSamHeader(string header_text) {
        auto hdr = sam_hdr_parse_(header_text.ptr, cast(int)header_text.length);
        if (cram_set_header(_fd, hdr) < 0)
            throw new CramException("Failed to set SAM header");
        if (_reference_fn !is null)
            cram_load_reference(_fd, _reference_fn.ptr);
        if (cram_write_SAM_hdr(_fd, _fd.header) < 0)
            throw new CramException("Failed to write SAM header");
    }

    private void convertToBamSeq(R)(ref R read, bam_seq_t* bam_seq) {
        with (bam_seq.core) {
            tid = read.ref_id;
            pos = read.position;
            bin = cast(ushort)read.bin.id;
            qual = read.mapping_quality;
            l_qname = cast(ubyte)(read.name.length + 1);
            flag = read.flag;
            n_cigar = cast(ushort)read.cigar.length;
            l_qseq = cast(int)read.sequence.length;
            mtid = read.mate_ref_id;
            mpos = read.mate_position;
            isize = read.template_length;
        }
        bam_seq.l_data = bam_seq.m_data = cast(int)(read.raw_data.length - 32);
        bam_seq.data = read.raw_data.ptr + 32;
    }

    void writeRecord(R)(auto ref R read) if(isBamRead!R) {
        // somewhat messy code because of the decision
        // to use ending zero byte of qname for flags
        auto offset = 32 + read.name.length;
        ubyte old_byte = read.raw_data[offset];
        read.raw_data[offset] = 0;
        scope(exit) read.raw_data[offset] = old_byte;

        bam_seq_t bam_seq;
        convertToBamSeq(read, &bam_seq);
        cram_put_bam_seq(_fd, &bam_seq);
    }

    void flush() {
        if (cram_flush(_fd) == -1)
            throw new CramException("Failed to flush block (" ~ _fn ~ ")");
    }

    void finish() {
        if (cram_close(_fd) == -1)
            throw new CramException("Failed to close CRAM file: " ~ _fn);
    }
}
