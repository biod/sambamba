module pileuprange;

import alignment;
import std.container;
import std.algorithm;

/// Represents a single pileup column
struct PileupColumn(R) {
    private {
        ulong _position;
        uint _coverage;
        R _reads;
    }

    this(ulong position, uint coverage, R reads) {
        _position = position;
        _coverage = coverage;
        _reads = reads;
    }

    /// Coverage at this position
    uint coverage() const @property {
        return _coverage;
    }

    /// Position on the reference
    ulong position() const @property {
        return _position;
    }

    /// Reads overlapping the position
    auto reads() @property {
        return _reads;
    }
}

auto pileupColumn(R)(ulong position, uint coverage, R reads) {
    return PileupColumn!R(position, coverage, reads);
}

/**
 * The class for iterating reference bases together with reads overlapping them.
 */
class PileupRange(R) {
    private {
        R _reads;
        uint _cur_pos;
        Tree _tree;

        static struct PostProcessedRead {
            Alignment read;
            uint end_position;
        }

        alias RedBlackTree!(PostProcessedRead, "a.end_position < b.end_position", true) Tree;
    }

    /**
     * Create new pileup iterator from a range of reads.
     */
    this(R reads) {
        _reads = reads;
        _tree = new Tree();
        if (!_reads.empty) {
            auto read = _reads.front;
            _cur_pos = read.position;
            _tree.stableInsert(PostProcessedRead(read, _cur_pos + read.basesCovered()));
            _reads.popFront();
            while (!_reads.empty) {
                read = _reads.front;
                if (read.position == _cur_pos) {
                    _tree.stableInsert(PostProcessedRead(read, _cur_pos + read.basesCovered()));
                    _reads.popFront();
                } else {
                    break;
                }
            }
        }
    }

    /// Returns PileupColumn struct corresponding to the current position.
    auto front() @property {
        return pileupColumn(cast(ulong)_cur_pos, cast(uint)_tree.length, map!"a.read"(_tree[]));
    }

    /// Whether all reads have been processed.
    auto empty() @property {
        return _reads.empty && _tree.empty;
    }

    /// Move to next position on the reference.
    void popFront() {
        ++_cur_pos;
        Alignment dummy;
        _tree.remove(_tree.lowerBound(PostProcessedRead(dummy, _cur_pos)));

        while (!_reads.empty && _reads.front.position == _cur_pos) {
            auto read = _reads.front;
            _tree.stableInsert(PostProcessedRead(read, read.position + read.basesCovered()));
            _reads.popFront();
        }
    }
}

/// Creates a pileup range from a range of reads.
auto pileup(R)(R reads) {
    return new PileupRange!R(reads);
}
