module rangetransformer;

import std.parallelism;
import std.traits : ReturnType;
import std.range : ElementType;

class RangeTransformer(alias func, R) {

    alias ElementType!R ElemR;

	alias ReturnType!func RetElem;

    alias Task!(func, ElemR) ETask;
    
    R range;
    
    uint num_of_threads;

    uint read = 0;
    uint consumed = 0;
   
    TaskPool pool;
    ETask[] buffer;

    private void initialize() {
        /*
         * Fills the buffer with first num_of_threads tasks
         * and puts them to pool so that they get executed
         * as we iterate we range.
         */
        foreach (i; 0..num_of_threads) {
            if (range.empty()) {
                break;
            }

            buffer[read] = scopedTask!func(range.front);
            pool.put(buffer[read]);
            range.popFront();

            ++read;
        }
    }

    this(R range, uint num_of_threads) {
        this.num_of_threads = num_of_threads;

        this.pool = new TaskPool(num_of_threads);

        this.range = range;
        
        buffer = new ETask[num_of_threads];

        initialize();
    }
   
    void finish() {
        pool.finish();
    }

    bool empty() @property {
        bool result = range.empty && read == consumed;
        if (result) {
            pool.finish();
        }
        return result;
    }

    RetElem front() @property {
        return buffer[consumed % num_of_threads].yieldForce();
    }

    void popFront() {
        if (range.empty) {
            ++consumed;
            return;
        }
        buffer[consumed % num_of_threads] = scopedTask!func(range.front);
        pool.put(buffer[consumed % num_of_threads]);
        ++read;
        range.popFront();
        ++consumed;
    }
}
