/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

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
/// $(P This module is used for iterating over columns of alignment.)
/// $(P The function makePileup is called on
/// a range of coordinate-sorted reads mapped to the same reference.
/// It returns an input range of columns.)
/// $(P This returned range can then be iterated with $(D foreach).
/// First column is located at the same position on the reference,
/// as the first base of the first read.
/// $(BR)
/// Each $(D popFront) operation advances current position on the
/// reference. The default behaviour is to exclude sites with zero coverage
/// from the iteration.)
/// $(P Each column keeps set of reads that overlap corresponding position
/// on the reference.
/// If reads contain MD tags, and makePileup was asked
/// to use them, reference base at the column is also available.)
/// $(BR)
/// Each read preserves all standard read properties
/// but also keeps column-related information, namely
/// <ul>
///     $(LI number of bases consumed from the read sequence so far)
///     $(LI current CIGAR operation and offset in it)
///     $(LI all CIGAR operations before and after current one)</ul>
/// $(BR)
/// It is clear from the above that current CIGAR operation cannot be an insertion.
/// The following are suggested ways to check for them:
/// <ul>
///     $(LI $(D cigar_after.length > 0 &&
///              cigar_operation_offset == cigar_operation.length - 1 &&
///              cigar_after[0].type == 'I'))
///     $(LI $(D cigar_before.length > 0 &&
///              cigar_operation_offset == 0 &&
///              cigar_before[$ - 1].type == 'I'))</ul>
/// $(BR)
/// Example:
/// ---------------------------------------------------------
/// import bio.std.hts.bam.reader, bio.std.hts.bam.pileup, std.stdio, std.algorithm : count;
/// void main() {
///     auto bam = new BamReader("file.bam");       // assume single reference and MD tags
///     auto pileup = bam.reads().makePileup(useMD);
///     foreach (column; pileup) {
///         auto matches = column.bases.count(column.reference_base);
///         if (matches < column.coverage * 2 / 3)
///             writeln(column.position);           // print positions of possible mismatches
///     }
/// }
/// ---------------------------------------------------------
module bio.std.hts.bam.pileup;

import bio.std.hts.bam.cigar;
import bio.std.hts.bam.read;
import bio.std.hts.bam.md.reconstruct;
import bio.std.hts.bam.splitter;

import std.algorithm;
import std.range;
import std.random;
import std.traits;
import std.conv;
import std.array;
import std.exception;

/// Represents a read aligned to a column
struct PileupRead(Read=bio.std.hts.bam.read.EagerBamRead) {

    Read read; ///
    alias read this;
    private alias read _read;

    /// Current CIGAR operation. One of 'M', '=', 'X', 'D', 'N.
    /// Use $(D cigar_after)/$(D cigar_before) to determine insertions.
    bio.std.hts.bam.read.CigarOperation cigar_operation() @property const {
        return _cur_op;
    }

    /// Number of bases consumed from the current CIGAR operation.
    uint cigar_operation_offset() @property const {
        return _cur_op_offset;
    }

    /// CIGAR operations after the current operation
    const(bio.std.hts.bam.read.CigarOperation)[] cigar_after() @property const {
        return _read.cigar[_cur_op_index + 1 .. $];
    }

    /// CIGAR operations before the current operation
    const(bio.std.hts.bam.read.CigarOperation)[] cigar_before() @property const {
        return _read.cigar[0 .. _cur_op_index];
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns read base
    /// at the current column. Otherwise, returns '-'.
    char current_base() @property const {
        assert(_query_offset <= _read.sequence_length);
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return _read.sequence[_query_offset];
        } else {
            return '-';
        }
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns
    /// Phred-scaled read base quality at the current column.
    /// Otherwise, returns 255.
    ubyte current_base_quality() @property const {
        assert(_query_offset <= _read.sequence_length);
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return _read.base_qualities[_query_offset];
        } else {
            return 255;
        }
    }

    /// Returns number of bases consumed from the read sequence.
    /// $(BR)
    /// More concisely,
    /// $(UL
    ///     $(LI if current CIGAR operation is 'M', '=', or 'X',
    ///       index of current read base in the read sequence)
    ///     $(LI if current CIGAR operation is 'D' or 'N',
    ///       index of read base after the deletion)
    /// )
    /// (in both cases indices are 0-based)
    int query_offset() @property const {
        assert(_query_offset <= _read.sequence_length);
        return _query_offset;
    }

    /// Returns duplicate
    PileupRead dup() @property {
        PileupRead r = void;
        r._read = _read; // logically const, thus no .dup here
        r._cur_op_index = _cur_op_index;
        r._cur_op = _cur_op;
        r._cur_op_offset = _cur_op_offset;
        r._query_offset = _query_offset;
        return r;
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
            _read = read;

            // find first M/=/X/D operation
            auto cigar = _read.cigar;
            for (_cur_op_index = 0; _cur_op_index < cigar.length; ++_cur_op_index) {
                _cur_op = cigar[_cur_op_index];
                if (_cur_op.is_reference_consuming) {
                    if (_cur_op.type != 'N') {
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

                auto cigar = _read.cigar;
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
            assert(_cur_op_index < _read.cigar.length, "Invalid read " ~ _read.name
                                                       ~ " - CIGAR " ~ _read.cigarString()
                                                       ~ ", sequence " ~ to!string(_read.sequence));
        }
    }
}

static assert(isBamRead!(PileupRead!BamRead));
//static assert(isBamRead!(PileupRead!(EagerBamRead!BamRead)));

/// Represents a single pileup column
struct PileupColumn(R) {
    private {
        ulong _position;
        int _ref_id = -1;
        R _reads;
        size_t _n_starting_here;
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

    /// Reads overlapping the position, sorted by coordinate
    auto reads() @property {
        return assumeSorted!compareCoordinates(_reads[]);
    }

    /// Reads that have leftmost mapped position at this column
    auto reads_starting_here() @property {
        return _reads[$ - _n_starting_here .. $];
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
    alias Unqual!(ElementType!R) Raw;
    alias EagerBamRead!Raw Eager;
    alias PileupRead!Eager Read;
    alias Read[] ReadArray;
    alias TColumn!ReadArray Column;

    private {
        R _reads;
        Column _column;
        Appender!ReadArray _read_buf;
        bool _skip_zero_coverage;
    }

    protected {
        // This is extracted into a method not only to reduce duplication
        // (not so much of it), but to allow to override it!
        // For that reason it is not marked as final. Overhead of virtual
        // function is negligible compared to computations in EagerBamRead
        // constructor together with inserting new element into appender.
        void add(ref Raw read) {
            _read_buf.put(PileupRead!Eager(Eager(read)));
        }
    }

    /**
     * Create new pileup iterator from a range of reads.
     */
    this(R reads, bool skip_zero_coverage) {
        _reads = reads;
        _read_buf = appender!ReadArray();
        _skip_zero_coverage = skip_zero_coverage;

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
                if (survived < i)
                {
                    data[survived] = data[i];
                }
                ++survived;
            }
        }

        for (size_t i = 0; i < survived; ++i) {
            data[i].incrementPosition();
        }
                                      // unless range is empty, this value is
        _read_buf.shrinkTo(survived);

        _column._n_starting_here = 0; // updated either in initNewReference()
                                      // or in the loop below

        if (!_reads.empty) {
            if (_reads.front.ref_id != _column._ref_id &&
                survived == 0) // processed all reads aligned to the previous reference
            {
                initNewReference();
            } else {
                size_t n = 0;
                while (!_reads.empty &&
                        _reads.front.position == pos &&
                        _reads.front.ref_id == _column._ref_id)
                {
                    auto read = _reads.front;
                    add(read);
                    _reads.popFront();
                    ++n;
                }
                _column._n_starting_here = n;

                // handle option of skipping sites with zero coverage
                if (survived == 0 && n == 0 && _skip_zero_coverage) {
                    // the name might be misleading but it does the trick
                    initNewReference();
                }
            }
        }

        _column._reads = _read_buf.data;
    }

    protected void initNewReference() {
        auto read = _reads.front;

        _column._position = read.position;
        _column._ref_id = read.ref_id;
        uint n = 1;
        add(read);

        _reads.popFront();

        while (!_reads.empty) {
            read = _reads.front;
            if (read.ref_id == _column.ref_id &&
                read.position == _column._position)
            {
                add(read);
                ++n;
                _reads.popFront();
            } else {
                break;
            }
        }

        _column._n_starting_here = n;
        _column._reads = _read_buf.data;
    }
}

/// Abstract pileup structure. S is type of column range.
struct AbstractPileup(R, S) {
    private R reads_;
    R reads() @property {
        return reads_;
    }

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

struct TakeUntil(alias pred, Range, Sentinel) if (isInputRange!Range)
{
    private Range _input;
    private Sentinel _sentinel;
    bool _done;

    this(Range input, Sentinel sentinel) {
        _input = input; _sentinel = sentinel; _done = _input.empty || predSatisfied();
    }

    @property bool empty() { return _done; }
    @property auto ref front() { return _input.front; }
    private bool predSatisfied() { return startsWith!pred(_input, _sentinel); }
    void popFront() { _input.popFront(); _done = _input.empty || predSatisfied(); }
}

auto takeUntil(alias pred, Range, Sentinel)(Range range, Sentinel sentinel) {
    return TakeUntil!(pred, Range, Sentinel)(range, sentinel);
}

auto pileupInstance(alias P, R)(R reads, ulong start_from, ulong end_at, bool skip_zero_coverage) {
    auto rs = filter!"a.basesCovered() > 0"(reads);
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
    auto sameref_rs = takeUntil!"a.ref_id != b"(rs, ref_id);
    alias typeof(sameref_rs) ReadRange;
    PileupRange!ReadRange columns = new P!ReadRange(sameref_rs, skip_zero_coverage);
    while (!columns.empty) {
        auto c = columns.front;
        if (c.position < start_from) {
            columns.popFront();
        } else {
            break;
        }
    }
    auto chopped = takeUntil!"a.position >= b"(columns, end_at);
    return AbstractPileup!(R, typeof(chopped))(reads, chopped, start_from, end_at, ref_id);
}

auto pileupColumns(R)(R reads, bool use_md_tag=false, bool skip_zero_coverage=true) {
    auto rs = filter!"a.basesCovered() > 0"(reads);
    alias typeof(rs) ReadRange;
    PileupRange!ReadRange columns;
    if (use_md_tag) {
        columns = new PileupRangeUsingMdTag!ReadRange(rs, skip_zero_coverage);
    } else {
        columns = new PileupRange!ReadRange(rs, skip_zero_coverage);
    }
    return columns;
}

/// Tracks current reference base
final static class PileupRangeUsingMdTag(R) :
    PileupRange!(R, PileupColumn)
{
    // The code is similar to that in reconstruct.d but here we can't make
    // an assumption about any particular read having non-zero length on reference.

    // current chunk of reference
    private alias typeof(_column._reads[].front) Read;
    private ReturnType!(dna!Read) _chunk;

    // end position of the current chunk on reference (assuming half-open interval)
    private uint _chunk_end_position;

    // next read from which we will extract reference chunk
    //
    // More precisely,
    // _next_chunk_provider = argmax (read => read.end_position)
    //                 {reads overlapping current column}
    private Read _next_chunk_provider;

    private bool _has_next_chunk_provider = false;

    // coverage at the previous location
    private ulong _prev_coverage;

    // we also track current reference ID
    private int _curr_ref_id = -1;

    ///
    this(R reads, bool skip_zero_coverage) {
        super(reads, skip_zero_coverage);
    }

    alias Unqual!(ElementType!R) Raw;

    //  Checks length of the newly added read and tracks the read which
    //  end position on the reference is the largest.
    //
    //  When reconstructed reference chunk will become empty, next one will be
    //  constructed from that read. This algorithm allows to minimize the number
    //  of reads for which MD tag will be decoded.
    protected override void add(ref Raw read) {
        // the behaviour depends on whether a new contig starts here or not
        bool had_zero_coverage = _prev_coverage == 0;

        super.add(read);

        // get wrapped read
        auto _read = _read_buf.data.back;

        // if we've just moved to another reference sequence, do the setup
        if (_read.ref_id != _curr_ref_id) {
            _curr_ref_id = _read.ref_id;

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
        _prev_coverage = 0;
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

    ///
    override void popFront() {
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
///     use_md_tag =  if true, use MD tag together with CIGAR
///                   to recover reference bases
///
///     start_from =  position from which to start
///
///     end_at     =  position before which to stop
///
/// $(BR)
/// That is, the range of positions is half-open interval
/// $(BR)
/// [max(start_from, first mapped read start position),
/// $(BR)
///  min(end_at, last mapped end position among all reads))
///
///     skip_zero_coverage = if true, skip sites with zero coverage
///
auto makePileup(R)(R reads,
                   bool use_md_tag=false,
                   ulong start_from=0,
                   ulong end_at=ulong.max,
                   bool skip_zero_coverage=true)
{
    if (use_md_tag) {
        return pileupInstance!PileupRangeUsingMdTag(reads, start_from, end_at, skip_zero_coverage);
    } else {
        return pileupInstance!PileupRange(reads, start_from, end_at, skip_zero_coverage);
    }
}

/// Allows to express the intention clearer.
enum useMD = true;

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

// multiple sequence alignment:
//                                                            ***
// ATTATGGACATTGTTTCCGTTATCATCATCATCATCATCATCATCATTATCATC
//       GACATTGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCAT---C
//          ATTGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCATCACC
//            TGTTTCCGTTATCATCATCATCATCATCATCATCATCATCATCATCAT---CACCAC
//                TCCGTTATCATCATCATCATCATCATCATCATCATCATCATCAT---CACCACCACC
//                   GTTATCATCATCATCATCATCATCATCATCATCATCATCAT---CATCGTCACCCTG
//                            ATCATCATCATAATCATCATCATCATCATCAT---CATCGTCACCCTGTGTTGAG
//                                      TCATCATCATCGTCAC------------------CCTGTGTTGAGGACAGAAGTAATTTCCCTTTCTTGGCT
//                                               TCATCATCATCATCACCACCACCACCCTGTGTTGAGGACAGAAGTAATATCCCT
//                                                            ---CACCACCACCCTGTGTTGAGGACAGAAGTAATTTCCCTTTCTTGGCTGGTCACC
//   *         *         *         *         *         *            *         *        *         *
//  760       770       780       790       800       810          820       830      840       850

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

    BamRead[] reads = new BamRead[10];

    foreach (i; iota(10)) {
        reads[i] = BamRead(readnames[i], sequences[i], cigars[i]);
        reads[i].position = positions[i];
        reads[i].ref_id = 0;
        reads[i]["MD"] = md_tags[i];
    }

    auto first_read_position = reads.front.position;
    auto reference = to!string(dna(reads));

    import std.stdio;
    // stderr.writeln("Testing pileup (low-level aspects)...");

    auto pileup = makePileup(reads, true, 796, 849, false);
    auto pileup2 = makePileup(reads, true, 0, ulong.max, false);
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
        assert(column.reference_base == reference[cast(size_t)(column.position - first_read_position)]);

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
            case 810:
                // last read is not yet fetched by pileup engine
                assert(column.reads[column.coverage - 2].cigar_after.front.type == 'D');
                pileup.popFront();
                break;
            case 817:
                assert(column.reads[column.coverage - 2].cigar_before.back.type == 'I');
                pileup.popFront();
                break;
            case 821:
                assert(column.reads[column.coverage - 3].cigar_operation.type == 'D');
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

    reads = [BamRead("r1", "CCCACATAGAAAGCTTGCTGTTTCTCTGTGGGAAGTTTTAACTTAGGTCAGCTT",
                       [CigarOperation(54, 'M')]),
             BamRead("r2", "TAGAAAGCTTGCTGTTTCTCTGTGGGAAGTTTTAACTTAGGTTAGCTTCATCTA",
                       [CigarOperation(54, 'M')]),
             BamRead("r3", "TTTTTCTTTCTTTCTTTGAAGAAGGCAGATTCCTGGTCCTGCCACTCAAATTTT",
                       [CigarOperation(54, 'M')]),
             BamRead("r4", "TTTCTTTCTTTCTTTGAAGAAGGCAGATTCCTGGTCCTGCCACTCAAATTTTCA",
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
                 map!(c => c.reference_base)(makePileup(reads, true, 0, ulong.max, false))));
}

struct PileupChunkRange(C) {
    private C _chunks;
    private ElementType!C _prev_chunk;
    private ElementType!C _current_chunk;
    private bool _empty;
    private ulong _beg = 0;
    private bool _use_md_tag;
    private ulong _start_from;
    private ulong _end_at;
    private int _chunk_right_end;

    private int computeRightEnd(ref ElementType!C chunk) {
        return chunk.map!(r => r.position + r.basesCovered()).reduce!max;
    }

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

                if (_current_chunk[0].ref_id < 0) continue;

                _beg = _current_chunk[0].position;
                if (_beg >= end_at) {
                    _empty = true;
                    break;
                }

                _chunk_right_end = computeRightEnd(_current_chunk);
                if (_chunk_right_end > start_from)
                    break;
            }
        }
    }

    bool empty() @property {
        return _empty;
    }

    auto front() @property {
        auto end_pos = _current_chunk[$-1].position;
        if (_chunks.empty || _chunks.front[0].ref_id != _current_chunk[$-1].ref_id)
            end_pos = _chunk_right_end;

        return makePileup(chain(_prev_chunk, _current_chunk),
                          _use_md_tag,
                          max(_beg, _start_from), min(end_pos, _end_at));
    }

    void popFront() {
        _prev_chunk = _current_chunk;

        while (true) {
            if (_chunks.empty) {
                _empty = true;
                return;
            }
            _current_chunk = _chunks.front;
            _chunks.popFront();

            if (_current_chunk[0].ref_id >= 0) break;
        }

        _chunk_right_end = computeRightEnd(_current_chunk);

        // if we changed reference, nullify prev_chunk
        if (_prev_chunk.length > 0 &&
            _prev_chunk[$ - 1].ref_id == _current_chunk[0].ref_id)
        {
            _beg = _prev_chunk[$-1].position;
        } else {
            _beg = _current_chunk[0].position;
            _prev_chunk.length = 0;
        }

        // keep only those reads in _prev_chunk that have overlap with the last one

        // 1) estimate read length
        enum sampleSize = 15;
        int[sampleSize] buf = void;
        int read_length = void;
        if (_prev_chunk.length <= sampleSize) {
            for (size_t k = 0; k < _prev_chunk.length; ++k) {
                buf[k] = _prev_chunk[k].sequence_length;
            }
            topN(buf[0.._prev_chunk.length], _prev_chunk.length / 2);
            read_length = buf[_prev_chunk.length / 2];
        } else {
            size_t i = 0;
            foreach (read; randomSample(_prev_chunk, sampleSize))
                buf[i++] = read.sequence_length;
            topN(buf[], sampleSize / 2);
            read_length = buf[sampleSize / 2];
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
            auto m = cast(size_t)(i + j) / 2;
            assert(m < _prev_chunk.length);
            auto p = _prev_chunk[m].position;
            if (p >= pos) {
                j = m - 1;
            } else {
                i = m + 1;
            }
        }

        _prev_chunk = _prev_chunk[cast(size_t)i .. $];
    }
}

/// This function constructs range of non-overlapping consecutive pileups from a range of reads
/// so that these pileups can be processed in parallel.
///
/// It's allowed to pass ranges of sorted reads with different ref. IDs,
/// they won't get mixed in any chunk.
///
/// Params:
///   use_md_tag =   recover reference bases from MD tag and CIGAR
///
///   block_size =   approximate amount of memory that each pileup will consume,
///                  given in bytes. (Usually consumption will be a bit higher.)
///
///   start_from =   position of the first column of the first chunk
///
///   end_at     =   position after the last column of the last chunk
///
/// $(BR)
/// WARNING:     block size should be big enough so that every block will share
///              some reads only with adjacent blocks.
///              $(BR)
///              As such, it is not recommended to reduce the $(I block_size).
///              But there might be a need to increase it in case of very high coverage.
auto pileupChunks(R)(R reads, bool use_md_tag=false, size_t block_size=16_384_000,
                     ulong start_from=0, ulong end_at=ulong.max) {
    auto chunks = chunksConsumingLessThan(reads, block_size);
    return PileupChunkRange!(typeof(chunks))(chunks, use_md_tag, start_from, end_at);
}
