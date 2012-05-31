import std.traits;
import std.typecons;
import std.parallelism;

class CacheDefaultImpl(R, Args...) {
    private R[Tuple!Args] cache; 

    R* lookup(Args args) {
        return tuple(args) in cache;
    }

    void associate(Args args, R value) {
        cache[tuple(args)] = value;
    }

    void update(Args args) {}
}

version(unittest) {
    /// keeps number of function evaluations
    static shared evaluations = 0;
}

auto memoize(alias func, 
             alias CacheImpl=CacheDefaultImpl, Args...)(Args args)
    if(isCallable!func && is(Args == ParameterTypeTuple!func)
                       && Args.length > 0) 
{
    alias ReturnType!func R; 

    static shared(CacheImpl!(R, Args)) cache;
    static shared bool init = false;

    if (!init) {
        synchronized {
            if (cache is null) {
                cache = new shared(CacheImpl!(R, Args))();
            }
            init = true;
        }
    }

    R* ret = (cast()cache).lookup(args);
    if (ret !is null) {
        return *ret;
    } else {
        synchronized(cache) {

            ret = (cast()cache).lookup(args);
            if (ret !is null) {
                return *ret;
            }

            debug {
                import std.stdio;
                writeln("evaluating... args: ", tuple(args));
            }

            version(unittest) {
                evaluations += 1;
            }

            auto result = func(args);
            (cast()cache).associate(args, result);
            (cast()cache).update(args);
            return result;
        }
    }
}

unittest {

    import core.thread;

    /// very simple function for testing
    int func(int x, int y) {
        return x * y;
    }

    /// 4 different argument values in total
    void thread_func() {
        memoize!func(5, 10);
        memoize!func(3, 7);
        memoize!func(6, 4);
        memoize!func(5, 10);
        memoize!func(7, 9);
    }
    
    Thread[] threads = new Thread[5];
    foreach (i; 0 .. 5) threads[i] = new Thread(&thread_func);
    foreach (i; 0 .. 5) threads[i].start();
    foreach (i; 0 .. 5) threads[i].join();

    /// number of evaluations must be the same as number of 
    /// different argument values (4 in this case)
    assert(evaluations == 4);
}

void main(){}
