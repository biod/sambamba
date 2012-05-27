/**
  Module for reading BAI files
*/
module bai.read;

import bai.chunk;
import bai.bin;

import std.stream;
import std.system;
import std.exception;
import std.algorithm;

struct Index {
    Bin[] bins;
    ulong[] ioffsets; // virtual file offsets of first alignments in intervals
}

struct BaiFile {
    Index[] indices;

    /// Initialize from stream which contains BAI data
    this(ref Stream stream) {
        _stream = stream;
        parse();
    }

    /// Open BAI file given either filename of BAM file or that of BAI file.
    this(string filename) {
        Stream fstream;
        if (endsWith(filename, ".bam")) {
            filename ~=  ".bai";
        }
        fstream = new BufferedFile(filename);
        Stream estream = new EndianStream(fstream, Endian.littleEndian);
        this(estream);
    }

private:
    Stream _stream;

    /// according to section 4.2 of SAM/BAM specification
    void parse() {
        auto magic = _stream.readString(4);
        enforce(magic == "BAI\1", "Invalid file format: expected BAI\\1");

        int n_ref;
        _stream.read(n_ref);
        indices.length = n_ref;

        foreach (i; 0 .. n_ref) {
            int n_bin;
            _stream.read(n_bin);
            indices[i].bins.length = n_bin;
            
            foreach (j; 0 .. n_bin) {
                _stream.read(indices[i].bins[j].id);

                int n_chunk;
                _stream.read(n_chunk);
                indices[i].bins[j].chunks.length = n_chunk;
                
                foreach (k; 0 .. n_chunk) {
                    _stream.read(indices[i].bins[j].chunks[k].beg);
                    _stream.read(indices[i].bins[j].chunks[k].end);
                }
            }

            int n_intv;
            _stream.read(n_intv);
            indices[i].ioffsets.length = n_intv;

            foreach (j; 0 .. n_intv) {
                _stream.read(indices[i].ioffsets[j]);
            }
        }
    }
}
