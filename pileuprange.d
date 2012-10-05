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
module pileuprange;

import alignment;
import reconstruct;

import std.algorithm;
import std.range;
import std.conv;
import std.array;
import std.exception;

static if (__traits(compiles, {import dcollections.LinkList;})) {
    import dcollections.LinkList;
} else {
    pragma(msg, "
    Pileup module uses containers from `dcollections` library.

    In order to install dcollections, execute following commands
    from the root directory of Sambamba:

    $ svn co http://svn.dsource.org/projects/dcollections/branches/d2
    $ cd d2
    $ ./build-lib-linux.sh
    $ cd ..
    $ mv d2/dcollections .
    $ mv d2/libdcollections.a .
    $ rm -rf d2/

    ");
}


/// Represents a read aligned to a column
struct PileupRead(Read=Alignment) {

    Read read;
    alias read this;
  
    /// Current CIGAR operation. Cannot be 'P' as padding operators are skipped.
    CigarOperation cigar_operation() @property const {
        return _cur_op;
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns read base
    /// at the current column. Otherwise, returns '-'.
    char current_base() @property const {
        assert(_query_offset <= read.sequence_length);
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return read.sequence[_query_offset];
        } else {
            return '-';
        }
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns
    /// Phred-scaled read base quality at the correct column.
    /// Otherwise, returns 255.
    ubyte current_base_quality() @property const {
        assert(_query_offset <= read.sequence_length);
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return read.phred_base_quality[_query_offset];
        } else {
            return 255;
        }
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns
    /// index of current base in the read sequence. 
    /// Otherwise, returns -1.
    int read_sequence_offset() @property const {
        assert(_query_offset <= read.sequence_length);
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return _query_offset;
        } else {
            return -1;
        }
    }

    private {    
        // index of current CIGAR operation in _read.cigar
        uint _cur_op_index;

        // current CIGAR operation
        CigarOperation _cur_op;

        // number of bases consumed from the current CIGAR operation
        uint _cur_op_offset;

        // number of bases consumed from the read sequence
        uint _query_offset;

        this(Read read) {
            this.read = read;

            // find first M/=/X/D operation
            auto cigar = read.cigar;
            for (_cur_op_index = 0; _cur_op_index < cigar.length; ++_cur_op_index) {
                assertCigarIndexIsValid();

                _cur_op = cigar[_cur_op_index];
                if (_cur_op.is_reference_consuming) {
                    if (_cur_op.operation != 'N') {     
                        break;
                    }
                } else if (_cur_op.is_query_consuming) {
                    _query_offset += _cur_op.length; // skip S and I operations
                }
            }

            assertCigarIndexIsValid();
        }

        // move one base to the right on the reference
        void incrementPosition() {
            ++_cur_op_offset;

            // if current CIGAR operation is D or N, query offset is untouched
            if (_cur_op.is_query_consuming) {
                ++_query_offset;
            }

            if (_cur_op_offset >= _cur_op.length) {

                _cur_op_offset = 0; // reset CIGAR operation offset

                auto cigar = read.cigar;
                // get next reference-consuming CIGAR operation (M/=/X/D/N)
                for (++_cur_op_index; _cur_op_index < cigar.length; ++_cur_op_index) {
                    _cur_op = cigar[_cur_op_index];
                    if (_cur_op.is_reference_consuming) {
                        break;
                    }

                    if (_cur_op.is_query_consuming) {
                        _query_offset += _cur_op.length;
                    }
                }

                assertCigarIndexIsValid();
            }
        }

        void assertCigarIndexIsValid() {
            assert(_cur_op_index < read.cigar.length, "Invalid read " ~ read.read_name 
                                                      ~ " - CIGAR " ~ read.cigarString()
                                                      ~ ", sequence " ~ to!string(read.sequence));
        }
    }
}

/// Represents a single pileup column
struct PileupColumn(R) {
    private {
        ulong _position;
        R _reads;
    }

    this(ulong position, R reads) {
        _position = position;
        _reads = reads;
    }

    /// Coverage at this position (equals to number of reads)
    ulong coverage() const @property {
        return _reads.length;
    }

    /// Position on the reference
    ulong position() const @property {
        return _position;
    }

    /// Reads overlapping the position
    auto reads() @property {
        return _reads[];
    }

    /// Shortcut for map!(read => read.current_base)(reads)
    auto bases() @property {
        return map!"a.current_base"(reads);
    }

    /// Shortcut for map!(read => read.current_base_quality)(reads)
    auto base_qualities() @property {
        return map!"a.current_base_quality"(reads);
    }

    /// Shortcut for map!(read => read.mapping_quality)(reads)
    auto read_qualities() @property {
        return map!"a.mapping_quality"(reads);
    }
}

/// Make a pileup column
auto pileupColumn(R)(ulong position, R reads) {
    return PileupColumn!R(position, reads);
}

/**
 * The class for iterating reference bases together with reads overlapping them.
 */
class PileupRange(R, alias TColumn=PileupColumn) {
    alias LinkList!(PileupRead!EagerAlignment) ReadList;
    alias TColumn!ReadList Column;

    private {
        R _reads;
        Column _column;
    }

    protected {
        // This is extracted into a method not only to reduce duplication
        // (not so much of it), but to allow to override it!
        // For that reason it is not marked as final. Overhead of virtual 
        // function is negligible compared to computations in EagerAlignment
        // constructor together with inserting a new node to the list.
        void add(ref Alignment read) {
            _column._reads.add(PileupRead!EagerAlignment(EagerAlignment(read)));
        }
    }

    /**
     * Create new pileup iterator from a range of reads.
     */
    this(R reads) {
        _reads = reads;

        _column._reads = new ReadList();

        if (!_reads.empty) {

            auto read = _reads.front;

            _column._position = read.position;
            add(read);

            _reads.popFront();

            while (!_reads.empty) {
                read = _reads.front;
                if (read.position == _column._position) {
                    add(read);
                    _reads.popFront();
                } else {
                    break;
                }
            }
        }
    }

    /// Returns PileupColumn struct corresponding to the current position.
    ref Column front() @property {
        return _column;
    }

    /// Whether all reads have been processed.
    bool empty() @property {
        return _reads.empty && _column._reads[].empty;
    }

    /// Move to next position on the reference.
    void popFront() {
        auto pos = ++_column._position;

        foreach (ref remove, ref _read; &_column._reads.purge) 
        {
            if (_read.end_position > pos) {
                _read.incrementPosition();
            } else {
                remove = true;
            }
        }
        
        while (!_reads.empty && _reads.front.position == pos) {
            auto read = _reads.front;
            add(read);
            _reads.popFront();
        }
    }
}

/// Creates a pileup range from a range of reads.
auto pileup(R)(R reads) {
    return new PileupRange!R(reads);
}

/// The same as pileup but allows to access bases of the reference.
/// That is, the current column has additional property reference_base().
///
/// NOTE: you can use this function only if reads have MD tags correctly set!
auto pileupWithReferenceBases(R)(R reads) {

    // The basic idea is:
    //   take PileupColumn, extend it using alias this (+reference_base());
    //   take PileupRange, extend it (override popFront()/add());
    //
    // What is the new functionality of add()? 
    //
    //  It will check length of the newly added read and track the read which
    //  end position on the reference is the largest. 
    //
    //  When reconstructed reference chunk will become empty, next one will be
    //  constructed from that read. This algorithm allows to minimize the number
    //  of reads for which MD tag will be decoded.
    //
    // What is the new functionality of popFront()?
    //
    //  It will pop an element off the reference chunk, and if it is empty,
    //  update the structure appropriately.

    static struct PileupColumnWithRefBases(R) {
        PileupColumn!R column;
        alias column this;

        this(ulong position, R reads) {
            column = PileupColumn!R(position, reads);
        }

        char reference_base() @property const {
            // zero coverage
            if (_reads[].empty) {
                return 'N';
            }

            // otherwise, the base must be set by PileupRangeWithRefBases code
            return _reference_base;
        }

        private char _reference_base;
    }

    final static class PileupRangeWithRefBases : 
        PileupRange!(R, PileupColumnWithRefBases) 
    {

        // The code is similar to that in reconstruct.d but here we can't make
        // an assumption about any particular read having non-zero length on reference.

        // current chunk of reference
        private typeof(dna(_reads.front)) _chunk;

        // end position of the current chunk on reference (assuming half-open interval)
        private uint _chunk_end_position;

        // next read from which we will extract reference chunk
        //
        // More precisely, 
        // _next_chunk_provider = argmax (read => read.end_position) 
        //                 {reads overlapping current column}
        private typeof(_column._reads[].front) _next_chunk_provider;

        private bool _has_next_chunk_provider = false;

        // coverage at the previous location
        private ulong _prev_coverage;

        this(R reads) {
            super(reads);

            if (!reads.empty) {
                auto _read = _column._reads[].back;

                // prepare first chunk
                _chunk = dna(_read.read.read); // two layers of wrapping 

                // set up _next_chunk_provider explicitly
                _next_chunk_provider = _read;
                _chunk_end_position = _read.end_position;
                _has_next_chunk_provider = true;

                _column._reference_base = _chunk.front;
                _chunk.popFront();
            }
        }

        protected override void add(ref Alignment read) {
            // the behaviour depends on whether a new contig starts here or not
            bool had_zero_coverage = _prev_coverage == 0;

            super.add(read);

            // get wrapped read
            auto _read = _column._reads[].back;

            // two subsequent next_chunk_providers must overlap
            // unless (!) there was a region with zero coverage in-between
            if (_read.position > _chunk_end_position && !had_zero_coverage) {
                return;
            }

            // compare with previous candidate and replace if this one is better
            if (_read.end_position > _chunk_end_position) {
                if (!_has_next_chunk_provider) {
                    _has_next_chunk_provider = true;
                    _next_chunk_provider = _read;
                } else if (_read.end_position > _next_chunk_provider.end_position) {
                    _next_chunk_provider = _read;
                }
            }
        }

        void popFront() {
            if (!_chunk.empty) {
                // update current reference base
                _column._reference_base = _chunk.front;

                _chunk.popFront();
            }

            // update _prev_coverage
            _prev_coverage = _column.coverage;

            // the order is important - maybe we will obtain new next_chunk_provider
            // during this call to popFront()
            super.popFront();

            // If we have consumed the whole current chunk,
            // we need to obtain the next one if it's possible.
            if (_chunk.empty && _has_next_chunk_provider) {
                _chunk = dna(_next_chunk_provider.read.read);

                debug {
                    
                    import std.stdio;
                    writeln();
                    writeln("position: ", _next_chunk_provider.position);
                    writeln("next chunk: ", to!string(_chunk));
                    
                }

                _chunk_end_position = _next_chunk_provider.end_position;

                _has_next_chunk_provider = false;

                _chunk.popFrontN(_column.position - _next_chunk_provider.position);

                _column._reference_base = _chunk.front;
                _chunk.popFront();
            }
        }
    }

    return new PileupRangeWithRefBases(reads);
}

unittest {
    import std.algorithm;
    import std.range;
    import std.array;

    // the set of reads below was taken from 1000 Genomes BAM file
    // NA20828.mapped.ILLUMINA.bwa.TSI.low_coverage.20101123.bam
    // (region 20:1127810-1127819)
    auto readnames = array(iota(10).map!(i => "r" ~ to!string(i))());

    auto sequences = ["ATTATGGACATTGTTTCCGTTATCATCATCATCATCATCATCATCATTATCATC",
                      "GACATTGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCATC",
                      "ATTGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCATCACC",
                      "TGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCATCACCAC",
                      "TCCGTTATCATCATCATCATCATCATCATCATCATCATCATCATCACCACCACC",
                      "GTTATCATCATCATCATCATCATCATCATCATCATCATCATCATCGTCACCCTG",
                      "TCATCATCATCATAATCATCATCATCATCATCATCATCGTCACCCTGTGTTGAG",
                      "TCATCATCATCGTCACCCTGTGTTGAGGACAGAAGTAATTTCCCTTTCTTGGCT",
                      "TCATCATCATCATCACCACCACCACCCTGTGTTGAGGACAGAAGTAATATCCCT",
                      "CACCACCACCCTGTGTTGAGGACAGAAGTAATTTCCCTTTCTTGGCTGGTCACC"];

    auto cigars = [[CigarOperation(54, 'M')],
                   [CigarOperation(54, 'M')],
                   [CigarOperation(50, 'M'), CigarOperation(3, 'I'), CigarOperation(1, 'M')],
                   [CigarOperation(54, 'M')],
                   [CigarOperation(54, 'M')],
                   [CigarOperation(54, 'M')],
                   [CigarOperation(2, 'S'), CigarOperation(52, 'M')],
                   [CigarOperation(16, 'M'), CigarOperation(15, 'D'), CigarOperation(38, 'M')],
                   [CigarOperation(13, 'M'), CigarOperation(3, 'I'), CigarOperation(38, 'M')],
                   [CigarOperation(54, 'M')]];

    auto positions = [758, 764, 767, 769, 773, 776, 785, 795, 804, 817];

    auto md_tags = ["47C6", "54", "51", "50T3", "46T7", "45A0C7", "11C24A0C14", 
                    "11A3T0^CATCATCATCACCAC38", "15T29T5", "2T45T5"];

    Alignment[] reads = array(iota(10).map!(i => Alignment(readnames[i], sequences[i], cigars[i]))());

    foreach (i; iota(10)) {
        reads[i].position = positions[i];
        reads[i].ref_id = 0;
        reads[i]["MD"] = md_tags[i];
    }

    auto first_read_position = reads.front.position;
    auto reference = to!string(dna(reads));

    import std.stdio;
    writeln("Testing pileup...");

    foreach (column; pileupWithReferenceBases(reads)) {
        // check that DNA is built correctly from MD tags and CIGAR
        assert(column.reference_base == reference[column.position - first_read_position]);

        switch (column.position) {
            case 796:
                assert(equal(column.bases, "CCCCCCAC"));
                break;
            case 805:
                assert(equal(column.bases, "TCCCCCCCC"));
                break;
            case 806:
                assert(equal(column.bases, "AAAAAAAGA"));
                break;
            case 821:
                assert(equal(column.bases, "AAGG-AA"));
                break;
            case 826:
                assert(equal(column.bases, "CCCCCC"));
                break;
            case 849:
                assert(equal(column.bases, "TAT"));
                break;
            default:
                break;
        }
    }

    // another set of reads, the same file, region 20:12360032-12360050
    // test the case when reference has some places with zero coverage

    reads = [Alignment("r1", "CCCACATAGAAAGCTTGCTGTTTCTCTGTGGGAAGTTTTAACTTAGGTCAGCTT",
                       [CigarOperation(54, 'M')]),
             Alignment("r2", "TAGAAAGCTTGCTGTTTCTCTGTGGGAAGTTTTAACTTAGGTTAGCTTCATCTA",
                       [CigarOperation(54, 'M')]),
             Alignment("r3", "TTTTTCTTTCTTTCTTTGAAGAAGGCAGATTCCTGGTCCTGCCACTCAAATTTT",
                       [CigarOperation(54, 'M')]),
             Alignment("r4", "TTTCTTTCTTTCTTTGAAGAAGGCAGATTCCTGGTCCTGCCACTCAAATTTTCA",
                       [CigarOperation(54, 'M')])];

    reads[0].position = 979;
    reads[0]["MD"] = "54";

    reads[1].position = 985;
    reads[1]["MD"] = "42C7C3";

    reads[2].position = 1046;
    reads[2]["MD"] = "54";

    reads[3].position = 1048;
    reads[3]["MD"] = "54";

    assert(equal(dna(reads), 
                 map!(c => c.reference_base)(pileupWithReferenceBases(reads))));
}
