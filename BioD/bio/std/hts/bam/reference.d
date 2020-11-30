/*
    This file is part of BioD.
    Copyright (C) 2012-2014    Artem Tarasov <lomereiter@gmail.com>

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
/// $(P Each BAM file contains reads aligned to different reference sequences.)
/// $(P These sequences have unique identifiers in BAM file, starting from 0.
/// Unmapped reads are associated with id = -1.)
/// $(P If BAI file is available, fast region queries are available, that is,
/// getting all reads that overlap given region. This is achieved via $(D opSlice) method.)
///
/// Example:
/// -----------------------------
/// import bio.std.hts.bam.reader, std.stdio;
/// ...
/// auto bam = new BamReader("file.bam");
/// auto refseq = bam["chr17"];
/// writeln(refseq.name, " - length ", refseq.length);
/// foreach (read; refseq[1234 .. 5678])
///     if (read.cigar.length > 1)
///         writeln(read.name, " ", read.cigarString());
/// -----------------------------
module bio.std.hts.bam.reference;

public import bio.std.hts.bam.referenceinfo;

import bio.std.hts.bam.readrange;
import bio.std.hts.bam.region;
import bio.std.hts.bam.randomaccessmanager;
import bio.core.bgzf.virtualoffset;

import contrib.undead.stream;
import std.exception;
import std.array;

///
struct ReferenceSequence {
    private int _ref_id;
   
    /// Name
    string name() @property const {
        return _info.name;
    }

    /// Length in base pairs
    int length() @property const {
        return _info.length;
    }

    /// Reference ID
    int id() @property const {
        return _ref_id;
    }

    /// Get alignments overlapping [start, end) region.
    /// $(BR)
    /// Coordinates are 0-based.
    auto opSlice(uint start, uint end) {
        enforce(start < end, "start must be less than end");
        enforce(_manager !is null, "random access is not available");
        enforce(_ref_id >= 0, "invalid reference id");
        return _manager.getReads(BamRegion(cast(uint)_ref_id, start, end));
    }

    /// Get all alignments for this reference
    auto opSlice() {
        return opSlice(0, length);
    }

    private alias typeof(opSlice().front) Read;
    private Read _first_read() @property {
        return opSlice().front.dup;
    }

    /// Virtual offset at which reads, aligned to this reference, start in BAM file.
    /// If there are no reads aligned to this reference, returns virtual
    /// offset of the EOF block if it's presented, or the end of file.
    bio.core.bgzf.virtualoffset.VirtualOffset startVirtualOffset() {
        auto reads = opSlice();
        if (reads.empty) {
            return _manager.eofVirtualOffset();
        }
        return reads.front.start_virtual_offset;
    }

    /// Virtual offset before which reads, aligned to this reference, stop.
    /// If there are no reads aligned to this reference, returns virtual
    /// offset of the EOF block if it's presented, or the end of file.
    bio.core.bgzf.virtualoffset.VirtualOffset endVirtualOffset() {

        if (opSlice().empty) {
            return _manager.eofVirtualOffset();
        }

        auto ioffsets = _manager.getBai().indices[_ref_id].ioffsets[];
        assert(ioffsets.length > 0);

        // Try to get startVirtualOffset of the next reference presented in the file.
        for (uint r = _ref_id + 1; r < _manager.getBai().indices.length; ++r) {
            auto reads = _manager.getReads(BamRegion(r, 0, uint.max));
            if (reads.empty) {
                continue;
            } else {
                return reads.front.start_virtual_offset;
            }
        }

        // However, this approach fails if there are unmapped reads coming after
        // this reference. We cannot just return _manager.eofVirtualOffset.

        auto last_offset = ioffsets[$ - 1];
        auto stream = _manager.createStreamStartingFrom(last_offset);
        auto last_few_reads = bamReadRange!withOffsets(stream, null);

        VirtualOffset result;
        assert(!last_few_reads.empty);
        foreach (read; last_few_reads) {
            if (read.ref_id == -1) break;
            result = read.end_virtual_offset;
        }

        return result;
    }
 
    /// First position on the reference overlapped by reads (0-based)
    /// $(BR)
    /// Returns -1 if set of reads is empty.
    int firstPosition() {
        auto reads = opSlice();
        if (reads.empty) {
            return -1;
        }
        return reads.front.position;
    }
   
    /// Last position on the reference overlapped by reads (0-based)
    /// $(BR)
    /// Returns -1 if set of reads is empty.
    int lastPosition() {
        // The key idea is
        //  1) use last offset from linear index
        //  2) loop through all remaining reads starting from there

        auto ioffsets = _manager.getBai().indices[_ref_id].ioffsets[];

        long index = ioffsets.length - 1;

        debug {
            int reads_processed = 0;
        }

        while (index >= 0) {
            auto offset = ioffsets[cast(size_t)index];

            auto stream = _manager.createStreamStartingFrom(offset);
            auto reads = bamReadRange(stream, null);

            int last_position = int.min;

            foreach (read; reads) {

                 debug {
                     reads_processed += 1;
                 }

                 if (read.ref_id != _ref_id) {
                     break;
                 }
                
                 if (read.position == -1) {
                     continue;
                 }

                 auto end_pos = read.position + read.basesCovered();
                 if (end_pos > last_position)
                     last_position = end_pos;
            }

            if (last_position != int.min) {
                debug {
                    import std.stdio;
                    stderr.writeln("[debug] ReferenceSequence.lastPosition() processed ",
                                   reads_processed, " reads");
                }
                return last_position - 1;
            }

            --index;
        }

        return firstPosition();
    }

    this(RandomAccessManager manager, int ref_id, ReferenceSequenceInfo info) {
        _manager = manager;
        _ref_id = ref_id;
        _info = info;
    }

private:
    RandomAccessManager _manager;
    ReferenceSequenceInfo _info;
}
