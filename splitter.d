module splitter;

import alignment;
import std.array;
import std.functional;
import std.range;
import std.traits;

/// Constructs range of chunks where sum { fn(read) | read in chunk }
/// does not exceed given number.
struct AlignmentRangeSplitter(R, alias fn) 
{
    this(R range, size_t threshold, bool split_by_ref) {
        _range = range;
        _size = threshold;
        _split_by_ref = split_by_ref;
        _appender = appender!(Alignment[])();
        getNextChunk();
    }

    private {
        R _range;
        bool _empty;
        bool _split_by_ref;
        size_t _size;
        Appender!(Alignment[]) _appender;
    }

    bool empty() @property {
        return _empty;
    }

    Alignment[] front() @property {
        return _appender.data.dup;
    }

    void popFront() {
        _appender.clear();
        getNextChunk(); 
    }

    private void getNextChunk() {
        if (_range.empty) {
            _empty = true;
            return;
        } 

        auto first_read = _range.front;
        _range.popFront();

        size_t total_size = first_read.size_in_bytes;
        auto average_size_estimate = unaryFun!fn(first_read);

        _appender.reserve(_size / average_size_estimate);
        _appender.put(first_read);

        while (total_size <= _size && !_range.empty) {
            auto read = _range.front;
            if (_split_by_ref && (read.ref_id != first_read.ref_id)) {
                break;
            }
            total_size += unaryFun!fn(read);
            _appender.put(read);
            _range.popFront();
        }
    }
}

/// Split range in chunks where total amount of memory consumed by all reads 
/// in the chunk is roughly chunk_size bytes.
///
/// Parameter $(D split_by_ref) specifies that each chunk should contain reads
/// aligned to the same reference. In most cases, this simplifies post-processing,
/// but in some cases this is not required, therefore it is optional.
auto chunksConsumingLessThan(R)(R alignments, size_t size_in_bytes, bool split_by_ref=true) {
    return AlignmentRangeSplitter!(R, "a.size_in_bytes")(alignments, size_in_bytes, split_by_ref);
}

/// Split range in chunks each containing no more than N reads
auto chunksOfSize(R)(R alignments, size_t N, bool split_by_ref=true) {
    return AlignmentRangeSplitter!(R, "1")(alignments, N, split_by_ref);
}
