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
import splitter;

import std.algorithm;
import std.range;
import std.random;
import std.conv;
import std.array;
import std.exception;

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
        int _ref_id = -1;
        R _reads;
    }

    /// Reference base. 'N' if not available.
    char reference_base() @property const {
        return _reference_base;
    }

    private char _reference_base = 'N';

    /// Coverage at this position (equals to number of reads)
    size_t coverage() const @property {
        return _reads.length;
    }

    /// Returns reference ID (-1 if unmapped)
    int ref_id() const @property {
        return _ref_id;
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

/**
 * The class for iterating reference bases together with reads overlapping them.
 */
class PileupRange(R, alias TColumn=PileupColumn) {
    alias PileupRead!EagerAlignment Read;
    alias Read[] ReadArray;
    alias TColumn!ReadArray Column;

    private {
        R _reads;
        Column _column;
        Appender!ReadArray _read_buf;
    }

    protected {
        // This is extracted into a method not only to reduce duplication
        // (not so much of it), but to allow to override it!
        // For that reason it is not marked as final. Overhead of virtual 
        // function is negligible compared to computations in EagerAlignment
        // constructor together with inserting a new node to the list.
        void add(ref Alignment read) {
            _read_buf.put(PileupRead!EagerAlignment(EagerAlignment(read)));
        }
    }

    /**
     * Create new pileup iterator from a range of reads.
     */
    this(R reads) {
        _reads = reads;
        _read_buf = appender!ReadArray();

        if (!_reads.empty) {
            initNewReference(); // C++ programmers, don't worry! Virtual tables in D
                                // are populated before constructor body is executed.
        }
    }

    /// Returns PileupColumn struct corresponding to the current position.
    ref Column front() @property {
        return _column;
    }

    /// Whether all reads have been processed.
    bool empty() @property {
        return _reads.empty && _read_buf.data.empty;
    }

    /// Move to next position on the reference.
    void popFront() {
        auto pos = ++_column._position;

        size_t survived = 0;
        auto data = _read_buf.data;
        for (size_t i = 0; i < data.length; ++i) {
            if (data[i].end_position > pos) {
                data[i].incrementPosition();
                data[survived++] = data[i];
            }
        }

        _read_buf.shrinkTo(survived);

        if (!_reads.empty) {
            if (_reads.front.ref_id != _column._ref_id && 
                survived == 0) // processed all reads aligned to the previous reference
            {
                initNewReference();
            } else {
                while (!_reads.empty && 
                        _reads.front.position == pos && 
                        _reads.front.ref_id == _column._ref_id) 
                {
                    auto read = _reads.front;
                    add(read);
                    _reads.popFront();
                }
            }
        }

        _column._reads = _read_buf.data;
    }

    protected void initNewReference() {
        auto read = _reads.front;

        _column._position = read.position;
        _column._ref_id = read.ref_id;
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

        _column._reads = _read_buf.data;
    }
}

/// Abstract pileup structure. S is type of column range.
struct AbstractPileup(S) {
    S columns;
    /// Pileup columns
    alias columns this;

    private {
        ulong _start_position;
        ulong _end_position;
        int _ref_id;
    }

    /// $(D start_from) parameter provided to a pileup function
    ulong start_position() @property const {
        return _start_position;
    }

    /// $(D end_at) parameter provided to a pileup function
    ulong end_position() @property const { 
        return _end_position;
    }

    /// Reference ID of all reads in this pileup.
    int ref_id() @property const {
        return _ref_id;
    }
}

auto pileupInstance(alias P, R)(R reads, ulong start_from, ulong end_at) {
    auto rs = filter!"!a.is_unmapped"(reads);
    while (!rs.empty) {
        auto r = rs.front;
        if (r.position + r.basesCovered() < start_from) {
            rs.popFront();
        } else {
            break;
        }
    }
    int ref_id = -1;
    if (!rs.empty) {
        ref_id = rs.front.ref_id;
    }
    auto sameref_rs = until!"a.ref_id != b"(rs, ref_id);
    alias typeof(sameref_rs) ReadRange;
    PileupRange!ReadRange columns = new P!ReadRange(sameref_rs);
    while (!columns.empty) {
        auto c = columns.front;
        if (c.position < start_from) {
            columns.popFront();
        } else {
            break;
        }
    }
    auto chopped = until!"a.position >= b"(columns, end_at);
    return AbstractPileup!(typeof(chopped))(chopped, start_from, end_at, ref_id);
}

auto pileupColumns(R)(R reads, bool use_md_tag=false) {
    auto rs = filter!"!a.is_unmapped"(reads);
    alias typeof(rs) ReadRange;
    PileupRange!ReadRange columns;
    if (use_md_tag) {
        columns = new PileupRangeUsingMdTag!ReadRange(rs);
    } else {
        columns = new PileupRange!ReadRange(rs);
    }
    return columns;
}

final static class PileupRangeUsingMdTag(R) : 
    PileupRange!(R, PileupColumn) 
{
    // The code is similar to that in reconstruct.d but here we can't make
    // an assumption about any particular read having non-zero length on reference.

    // current chunk of reference
    private typeof(dna(_column._reads[].front)) _chunk;

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

    // we also track current reference ID
    private int _curr_ref_id = -1;

    this(R reads) {
        super(reads);
    }

    //  Checks length of the newly added read and tracks the read which
    //  end position on the reference is the largest. 
    //
    //  When reconstructed reference chunk will become empty, next one will be
    //  constructed from that read. This algorithm allows to minimize the number
    //  of reads for which MD tag will be decoded.
    protected override void add(ref Alignment read) {
        // the behaviour depends on whether a new contig starts here or not
        bool had_zero_coverage = _prev_coverage == 0;

        super.add(read);

        // get wrapped read
        auto _read = _read_buf.data.back;

        // if we've just moved to another reference sequence, do the setup
        if (_read.ref_id != _curr_ref_id) {
            _curr_ref_id = _read.ref_id;

            _prev_coverage = 0;
            _has_next_chunk_provider = true;
            _next_chunk_provider = _read;
            return;
        }

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

    protected override void initNewReference() {
        super.initNewReference();
        if (_has_next_chunk_provider) {
            // prepare first chunk
            _chunk = dna(_next_chunk_provider);
            _chunk_end_position = _next_chunk_provider.end_position;
            _has_next_chunk_provider = false;
            _column._reference_base = _chunk.front;
            _chunk.popFront();
        } else {
            _column._reference_base = 'N';
        }
    }

    void popFront() {
        if (!_chunk.empty) {
            // update current reference base
            _column._reference_base = _chunk.front;

            _chunk.popFront();
        } else {
            _column._reference_base = 'N';
        }

        // update _prev_coverage
        _prev_coverage = _column.coverage;

        // the order is important - maybe we will obtain new next_chunk_provider
        // during this call to popFront()
        super.popFront();

        // If we have consumed the whole current chunk,
        // we need to obtain the next one if it's possible.
        if (_chunk.empty && _has_next_chunk_provider) {
            _chunk = dna(_next_chunk_provider);

            debug {
            /*  import std.stdio;
                writeln();
                writeln("position: ", _next_chunk_provider.position);
                writeln("next chunk: ", to!string(_chunk));
                */
            }

            _chunk_end_position = _next_chunk_provider.end_position;

            _has_next_chunk_provider = false;

            _chunk.popFrontN(cast(size_t)(_column.position - _next_chunk_provider.position));

            _column._reference_base = _chunk.front;
            _chunk.popFront();
        }
    }
}

/// Creates a pileup range from a range of reads.
/// Note that all reads must be aligned to the same reference.
///
/// See $(D PileupColumn) documentation for description of range elements.
/// Note also that you can't use $(D std.array.array()) function on pileup
/// because it won't make deep copies of underlying data structure.
/// (One might argue that in this case it would be better to use opApply,
/// but typically one would use $(D std.algorithm.map) on pileup columns
/// to obtain some numeric characteristics.)
/// 
/// Params:
///
///     use_md_tag -  if true, use MD tag together with CIGAR
///                   to recover reference bases
///
///     start_from -  position from which to start
///
///     end_at     -  position before which to stop
///
/// That is, the range of positions is half-open interval 
/// [max(start_from, first mapped read start position), 
///  min(end_at, last mapped read end position))
auto makePileup(R)(R reads, 
                   bool use_md_tag=false,
                   ulong start_from=0, 
                   ulong end_at=ulong.max) 
{
    if (use_md_tag) {
        return pileupInstance!PileupRangeUsingMdTag(reads, start_from, end_at);
    } else {
        return pileupInstance!PileupRange(reads, start_from, end_at);
    }
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

    Alignment[] reads = new Alignment[10];

    foreach (i; iota(10)) {
        reads[i] = Alignment(readnames[i], sequences[i], cigars[i]);
        reads[i].position = positions[i];
        reads[i].ref_id = 0;
        reads[i]["MD"] = md_tags[i];
    }

    auto first_read_position = reads.front.position;
    auto reference = to!string(dna(reads));

    import std.stdio;
    writeln("Testing pileup (low-level aspects)...");

    auto pileup = makePileup(reads, true, 796, 849);
    auto pileup2 = makePileup(reads, true);
    assert(pileup.front.position == 796);
    assert(pileup.start_position == 796);
    assert(pileup.end_position == 849);

    while (pileup2.front.position != 796) {
        pileup2.popFront();
    }

    while (!pileup.empty) {
        auto column = pileup.front;
        auto column2 = pileup2.front;
        assert(column.coverage == column2.coverage);
        pileup2.popFront();

        // check that DNA is built correctly from MD tags and CIGAR
        assert(column.reference_base == reference[column.position - first_read_position]);

        switch (column.position) {
            case 796:
                assert(equal(column.bases, "CCCCCCAC"));
                pileup.popFront();
                break;
            case 805:
                assert(equal(column.bases, "TCCCCCCCC"));
                pileup.popFront();
                break;
            case 806:
                assert(equal(column.bases, "AAAAAAAGA"));
                pileup.popFront();
                break;
            case 821:
                assert(equal(column.bases, "AAGG-AA"));
                pileup.popFront();
                break;
            case 826:
                assert(equal(column.bases, "CCCCCC"));
                pileup.popFront();
                break;
            case 849:
                assert(equal(column.bases, "TAT"));
                pileup.popFront();
                assert(pileup.empty);
                break;
            default:
                pileup.popFront();
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
    reads[0].ref_id = 0;

    reads[1].position = 985;
    reads[1]["MD"] = "42C7C3";
    reads[1].ref_id = 0;

    reads[2].position = 1046;
    reads[2]["MD"] = "54";
    reads[2].ref_id = 0;

    reads[3].position = 1048;
    reads[3]["MD"] = "54";
    reads[3].ref_id = 0;

    assert(equal(dna(reads), 
                 map!(c => c.reference_base)(makePileup(reads, true))));
}

/// However, it's not effective to do operations on pileup in a single thread.
/// This function constructs range of consecutive pileups from a range of reads
/// so that these pileups can be processed in parallel.
/// 
/// It's allowed to pass ranges of sorted reads with different ref. IDs,
/// they won't get mixed in any chunk.
/// 
/// Params:
///   use_md_tag -   recover reference bases from MD tag and CIGAR
///
///   block_size -   approximate amount of memory that each pileup will consume,
///                  given in bytes. (Usually consumption will be a bit higher.)
///
///   start_from -   position of the first column of the first chunk
///
///   end_at     -   position after the last column of the last chunk
///
/// @@WARNING@@: block size should be big enough so that every block will share
///              some reads only with adjacent blocks. 
///              As an example, you should not use block size 10000 with 454 reads,
///              since each block will contain only about ~10 alignment records.
///              The rule of thumb is 
///              block size = 4 * largest coverage * average read length
///              (e.g. for 454 reads, max. coverage 500 you would use size 1_000_000)
auto pileupChunks(R)(R reads, bool use_md_tag=false, size_t block_size=16_384_000,
                     ulong start_from=0, ulong end_at=ulong.max) {
    auto chunks = chunksConsumingLessThan(reads, block_size);

    static struct Result(C) {
        private C _chunks;
        private Alignment[] _prev_chunk;
        private Alignment[] _current_chunk;
        private bool _empty;
        private ulong _beg = 0;
        private bool _use_md_tag;
        private ulong _start_from;
        private ulong _end_at;

        this(C chunks, bool use_md_tag, ulong start_from, ulong end_at) {
            _chunks = chunks; 
            _use_md_tag = use_md_tag;
            _start_from = start_from;
            _end_at = end_at;
            while (true) {
                if (_chunks.empty) {
                    _empty = true;
                } else {
                    _current_chunk = _chunks.front;
                    _chunks.popFront();

                    if (_beg >= end_at) {
                        _empty = true;
                        break;
                    }

                    auto last_read = _current_chunk[$-1];
                    if (last_read.position + last_read.basesCovered() > start_from) {
                        break;
                    }
                }
            }
        }

        bool empty() @property {
            return _empty;
        }

        auto front() @property {
            return makePileup(chain(_prev_chunk, _current_chunk), 
                              _use_md_tag,
                              max(_beg, _start_from), 
                              min(_current_chunk[$-1].position, _end_at));
        }

        void popFront() {
            _prev_chunk = _current_chunk;

            if (_chunks.empty) {
                _empty = true;
                return;
            }
            _current_chunk = _chunks.front;
            _chunks.popFront();

            assert(_prev_chunk.length > 0);
            _beg = _prev_chunk[$-1].position;

            // keep only those reads in _prev_chunk that have overlap with the last one
            
            // 1) estimate read length
            int[15] buf = void;
            int read_length = void;
            if (_prev_chunk.length <= 15) {
                for (size_t k = 0; k < _prev_chunk.length; ++k) {
                    buf[k] = _prev_chunk[k].sequence_length;
                }
                topN(buf[0.._prev_chunk.length], _prev_chunk.length / 2);
                read_length = buf[_prev_chunk.length / 2];
            } else {
                copy(map!"a.sequence_length"(randomSample(_prev_chunk, 15)), buf[]);
                topN(buf[], 7);
                read_length = buf[7];
                debug {
                    import std.stdio;
                    stderr.writeln("[pileupChunks] read_length=", read_length);
                }
            }

            // 2) do binary search for those reads that start from (_beg - 2 * read_length)
            //    (it's an experimental fact that almost none of reads consumes that much
            //     on a reference sequence)
            auto pos = _beg - 2 * read_length;
            long i = 0;
            long j = _prev_chunk.length - 1;
            // positions of _prev_chunk[0 .. i] are less than pos,
            // positions of _prev_chunk[j + 1 .. $] are more or equal to pos.
            
            while (i <= j) { 
                auto m = (i + j) / 2;
                assert(m < _prev_chunk.length);
                auto p = _prev_chunk[m].position;
                if (p >= pos) {
                    j = m - 1;
                } else {
                    i = m + 1;
                }
            }

            _prev_chunk = _prev_chunk[i .. $];
        }
    }

    return Result!(typeof(chunks))(chunks, use_md_tag, start_from, end_at);
}
