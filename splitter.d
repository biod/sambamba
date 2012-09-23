module splitter;

import alignment;
import std.array;
import std.functional;
import std.range;
import std.traits;

/// Constructs range of chunks where sum { fn(read) | read in chunk }
/// does not exceed given number.
struct AlignmentRangeSplitter(R, alias fn) 
    if (isInputRange!R && is(Unqual!(ElementType!R) == Alignment))
{
    this(R range, size_t threshold) {
        _range = range;
        _size = threshold;
        _appender = appender!(Alignment[])();
        getNextChunk();
    }

    private {
        R _range;
        bool _empty;
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
            total_size += unaryFun!fn(read);
            _appender.put(read);
            _range.popFront();
        }
        debug {
            import std.stdio;
            stderr.writeln(total_size);
        }
    }
}

/// Split range in chunks where total amount of memory consumed by all reads 
/// in the chunk is roughly chunk_size bytes.
auto chunksConsumingLessThan(R)(R alignments, size_t size_in_bytes) {
    return AlignmentRangeSplitter!(R, "a.size_in_bytes")(alignments, size_in_bytes);
}

/// Split range in chunks each containing no more than N reads
auto chunksOfSize(R)(R alignments, size_t N) {
    return AlignmentRangeSplitter!(R, "1")(alignments, N);
}
