module utils.inputrangechunks;

import std.range;
import std.exception : enforce;

class Chunks(R)
    if(isInputRange!R)
{

private:
    R range;
    alias ElementType!R E;
    uint chunk_size;

    E[] buffer;

    void fillBuffer() {
        buffer = new E[chunk_size];
        for (auto i = 0; i < chunk_size; i++) {
            if (range.empty) {
                buffer.length = i;
                break;
            }
            buffer[i] = range.front;
            range.popFront();
        }
    }

public:

    this(R range, uint chunk_size) {
        enforce(chunk_size > 0);
        this.range = range;
        this.chunk_size = chunk_size; 
        fillBuffer();
    }

    bool empty() @property {
        return buffer.length == 0;
    }

    E[] front() @property {
        return buffer;    
    }

    void popFront() {
        fillBuffer();
    }
}

auto chunkedInputRange(R)(R range, uint chunk_size) {
    return new Chunks!R(range, chunk_size);
}

unittest {
    import std.algorithm;
    assert(equal(chunkedInputRange(iota(1, 6), 2), [[1, 2], [3, 4], [5]]));
    assert(equal(chunkedInputRange(iota(1, 7), 2), [[1, 2], [3, 4], [5, 6]]));
    assert(equal(chunkedInputRange([1], 10), [[1]]));
    auto r = iota(25);
    import std.stdio;
    writeln(joiner(chunkedInputRange(r, 7)));
    assert(equal(joiner(chunkedInputRange(r, 7)), r));
}
