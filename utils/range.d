module utils.range;

import std.range;
import std.exception;
import std.algorithm;
import std.parallelism;
import std.functional;
import std.array;

/// Keeps a cyclic buffer of size $(D amount)
/// which is filled at the construction.
/// After that, each popFront() is accompanied
/// by fetching next element from the original range.
///
/// The function is useful when, for instance, range of Tasks
/// is being decorated, because it allows to keep a certain amount
/// of them being executed simultaneously, utilizing all
/// CPU cores.
auto prefetch(Range)(Range r, uint amount) {

    enforce(amount > 0, "Amount of elements to prefetch must be positive");

    struct Result {
        alias ElementType!Range E;

        this(Range range, uint amount) {
            _elements = new E[amount];
            _range = range;
            _amount = amount;
            foreach (i; 0 .. _amount) {
                if (_range.empty) {
                    break;
                }
                _elements[i] = _range.front;
                ++ _read;
                _range.popFront();
            }
        }
        
        bool empty() @property {
            return _range.empty() && _read == _consumed;
        }
        
        E front() @property {
            return _elements[_consumed % _amount];
        }

        void popFront() @property {
            if (_range.empty) {
                ++ _consumed;
                return;
            }

            _elements[_consumed % _amount] = _range.front;
            ++ _consumed;

            _range.popFront();
            ++ _read;
        }
    private:
        Range _range;
        int _read = 0;
        int _consumed = 0;
        uint _amount;
        E[] _elements = void;
    }

    return Result(r, amount);
}

unittest {
    import std.algorithm;

    ubyte[] emptyrange = [];
    assert(equal(emptyrange, prefetch(emptyrange, 42)));

    auto range = [1, 2, 3, 4, 5];
    assert(equal(range, prefetch(range, 1)));
    assert(equal(range, prefetch(range, 3)));
    assert(equal(range, prefetch(range, 5)));
    assert(equal(range, prefetch(range, 7)));
}

/// Takes arbitrary input range as an input and returns
/// another range which produces arrays of original elements
/// of size $(D chunk_size).
///
/// Useful for setting granularity in parallel applications.
/// $(D std.algorithm.joiner) composed with $(D chunked) 
/// produces same elements as were in the original range.
/// 
/// The difference from $(D std.range.chunks) is that
/// any input range is allowed, no slicing or length is required.
/// The cost is memory allocations for chunks.
auto chunked(R)(R range, uint chunk_size) {

	struct Result {

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
	}

    return Result(range, chunk_size);
}

unittest {
    import std.algorithm;

    assert(equal(chunked(iota(1, 6), 2), [[1, 2], [3, 4], [5]]));
    assert(equal(chunked(iota(1, 7), 2), [[1, 2], [3, 4], [5, 6]]));
    assert(equal(chunked([1], 10), [[1]]));
	assert(equal(chunked(iota(1, 10), 7), [[1, 2, 3, 4, 5, 6, 7], [8,9]]));

    auto r = iota(25);
    assert(equal(joiner(chunked(r, 7)), r));
} 

/// Version of parallel map using cyclic buffer with prefetching.
/// Uses combination of chunked, prefetch, joiner, and std.parallelism
/// 
/// Params:
/// 	prefetch_amount -   how many chunks will be prefetched
///     chunk_size 		-   the maximum size of each chunk
auto parallelTransform(alias func, Range)(Range r, 
										  uint chunk_size=1, 
										  uint prefetch_amount=totalCPUs-1)
{
	alias ElementType!Range E;

	static auto createTask(E[] elements) {
		auto task =  task!(pipe!(map!(unaryFun!func), array))(elements);
		taskPool.put(task);
		return task;
	}

	auto chunks = chunked(r, chunk_size);
	auto tasks = map!createTask(chunks);
	auto prefetched = prefetch(tasks, prefetch_amount);
	return joiner(map!"a.yieldForce()"(prefetched));
}

unittest {
	auto range = iota(100);
	assert(equal(parallelTransform!"a * a"(range), map!"a * a"(range)));
} 
