module pileuprange;

import alignment;
import std.algorithm;
import std.conv;
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
    /// at the current column. Otherwise, returns 'N'.
    char base() @property const {
        if (_cur_op.is_query_consuming && _cur_op.is_reference_consuming) {
            return read.sequence[_query_offset];
        } else {
            return 'N';
        }
    }

    /// If current CIGAR operation is one of 'M', '=', or 'X', returns
    /// Phred-scaled read base quality at the correct column.
    /// Otherwise, returns 255.
    ubyte quality() @property const {
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
            ++_query_offset;

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
    auto reads() @property const {
        return _reads;
    }
}

/// Make a pileup column
auto pileupColumn(R)(ulong position, R reads) {
    return PileupColumn!R(position, reads);
}

/**
 * The class for iterating reference bases together with reads overlapping them.
 */
class PileupRange(R) {
    alias LinkList!(PileupRead!EagerAlignment) ReadList;
    alias PileupColumn!ReadList Column;

    private {
        R _reads;
        Column _column;
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
            _column._reads.add(PileupRead!EagerAlignment(EagerAlignment(read)));

            _reads.popFront();

            while (!_reads.empty) {
                read = _reads.front;
                if (read.position == _column._position) {
                    _column._reads.add(PileupRead!EagerAlignment(EagerAlignment(read)));

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
            _column._reads.add(PileupRead!EagerAlignment(EagerAlignment(read)));
            _reads.popFront();
        }
    }
}

/// Creates a pileup range from a range of reads.
auto pileup(R)(R reads) {
    return new PileupRange!R(reads);
}
