/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

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
module bio.std.hts.bam.baifile;

public import bio.core.bgzf.chunk;
public import bio.std.hts.bam.bai.bin;
import bio.core.bgzf.virtualoffset;
import bio.std.hts.bam.constants;

//TODO Remove stream dependency!
import contrib.undead.stream;
import std.system;
import std.exception;
import std.algorithm;
import std.conv;
import std.range;
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

    /// Get (approximate) virtual offset of the first alignment overlapping $(D position)
    /// 
    /// Returned virtual offset is less or equal to real offset.
    VirtualOffset getMinimumOffset(int position) {
        int pos = max(0, position);
        int _i = min(pos / BAI_LINEAR_INDEX_WINDOW_SIZE, cast(int)ioffsets.length - 1);
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

        if (!endsWith(filename, ".bai")) {
            /// Unfortunately, std.path.addExt is going to be deprecated

            auto first_filename = filename ~ ".bai";
            auto second_filename = to!string(retro(find(retro(filename), '.'))) ~ "bai";

            if (std.file.exists(first_filename)) {
                fstream = new BufferedFile(absolutePath(first_filename));
            } else {
                if (std.file.exists(second_filename)) {
                    fstream = new BufferedFile(absolutePath(second_filename));
                } else {
                    throw new Exception("searched for " ~ first_filename ~ " or " ~
                                        second_filename ~ ", found neither");
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
        indices = uninitializedArray!(Index[])(n_ref);

        foreach (i; 0 .. n_ref) {
            int n_bin = void;
            _stream.read(n_bin);
            indices[i].bins = uninitializedArray!(Bin[])(n_bin);
            
            foreach (j; 0 .. n_bin) {
                uint id = void;
                _stream.read(id);
                auto bin = Bin(id);

                int n_chunk = void;
                _stream.read(n_chunk);
                bin.chunks = uninitializedArray!(Chunk[])(n_chunk);
                
                foreach (k; 0 .. n_chunk) {
                    ulong tmp = void;
                    _stream.read(tmp);
                    bin.chunks[k].beg = VirtualOffset(tmp);
                    _stream.read(tmp);
                    bin.chunks[k].end = VirtualOffset(tmp);
                }
                
                indices[i].bins[j] = bin;
            }

            int n_intv = void;
            _stream.read(n_intv);
            indices[i].ioffsets = uninitializedArray!(VirtualOffset[])(n_intv);

            foreach (j; 0 .. n_intv) {
                ulong tmp = void;
                _stream.read(tmp);
                indices[i].ioffsets[j] = VirtualOffset(tmp);
            }
        }
    }
}
