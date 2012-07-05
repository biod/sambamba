/**
  Module for reading BAI files
*/
module bai.read;

import bai.chunk;
import bai.bin;
import virtualoffset;
import constants;

import std.stream;
import std.system;
import std.exception;
import std.algorithm;
import std.file;
import std.path;

/// Represents index for a single reference
struct Index {
    /// Information about bins
    Bin[] bins; 

    /// Virtual file offsets of first alignments overlapping 16384-byte windows
    /// on the reference sequence. This linear index is used to reduce amount
    /// of file seeks for region queries, since with its help one can reduce the
    /// number of chunks to be investigated based on their end position.
    ///
    ///
    /// Suppose you have a region [beg, end) and want to do a query.
    ///
    /// Here's the reference:
    /// [....................!..............!.................................]
    ///                     beg            end
    ///
    /// Here's the same reference with 16384-byte long windows:
    /// [%...........%.......!....%.........!..%...........%...........%......]
    ///                     beg            end
    /// [ 1st window][ 2nd window][...
    ///
    /// With linear index, we can take the second window, find out what is 
    /// the minimum virtual offset among alignments overlapping this window,
    /// and skip all chunks which end position is less or equal to this offset:
    ///
    /// [........@...........!..............!.................................]
    ///   .  ..min. offset   beg           end
    ///   [  ).        .                              <- this chunk is skipped
    ///       [        )                              <- this one is not
    ///
    VirtualOffset[] ioffsets; 

    // Get virtual offset of the first alignment overlapping $(D position)
    VirtualOffset getMinimumOffset(int position) {
        int pos = max(0, position);
        int _i = min(pos >> BAI_LINEAR_INDEX_SHIFT, cast(int)ioffsets.length - 1);
        auto min_offset = (_i == -1) ? VirtualOffset(0) : ioffsets[_i];
        return min_offset;
    }
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
            /// Unfortunately, std.path.addExt is going to be deprecated
            if (std.file.exists(filename ~ ".bai")) {
                fstream = new BufferedFile(absolutePath(filename ~ ".bai"));
            } else {
                filename = filename[0 .. $ - 3] ~ "bai";
                if (std.file.exists(filename)) {
                    fstream = new BufferedFile(absolutePath(filename));
                } else {
                    throw new Exception("searched for " ~ filename ~ " or " ~
                                        filename[0..$-1] ~ "m.bai" ~ ", found neither");
                }
            }
        } else {
            fstream = new BufferedFile(filename);
        }

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
                    ulong tmp;
                    _stream.read(tmp);
                    indices[i].bins[j].chunks[k].beg = VirtualOffset(tmp);
                    _stream.read(tmp);
                    indices[i].bins[j].chunks[k].end = VirtualOffset(tmp);
                }
            }

            int n_intv;
            _stream.read(n_intv);
            indices[i].ioffsets.length = n_intv;

            foreach (j; 0 .. n_intv) {
                ulong tmp;
                _stream.read(tmp);
                indices[i].ioffsets[j] = VirtualOffset(tmp);
            }
        }
    }
}
