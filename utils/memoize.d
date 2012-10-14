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
module utils.memoize;

import std.traits;
import std.typecons;

debug {
    import std.stdio;
}

/// Elements are not removed from the cache.
class BasicCache(uint maxItems, K, V) {
    private V[K] cache; 

    V* lookup(K key) {
        return key in cache;
    }

    void put(K key, V value) {
        cache[key] = value;
    }
}

/// First-in first-out element removal.
class FifoCache(uint maxItems, K, V) {
    private V[K] cache; 
    private int items = 0;
    bool full = false;
    private K[maxItems] keys; // cyclic queue
    private uint removals;

    V* lookup(K key) {
        return key in cache;
    }

    void put(K key, V value) {
        cache[key] = value;

        if (full) {
            cache.remove(keys[items]);
            removals += 1;

            if (removals % maxItems == 0) {
                cache.rehash;
            }
        }

        keys[items] = key;
        items += 1;

        if (items == maxItems) {
            full = true;
            items = 0;
        }
    }
}

/// Least used strategy
class LuCache(uint maxItems, K, V) {
    private {
        V[K] cache;
        int[K] counter;
        uint removals;
    }

    V* lookup(K key) {
        auto result = key in cache;
        if (result !is null) {
            counter[key] += 1;
        }
        return result;
    }

    void put(K key, V value) {
        if (counter.length >= maxItems) {
            // delete one element before inserting next
            int min = int.max;
            K min_key;
            foreach (k, i; counter) {
                if (i < min) {
                    min = i;
                    min_key = k;
                }
            }
            
            cache.remove(min_key);
            counter.remove(min_key);
            removals += 1;

            if (removals % maxItems == 0) {
                cache.rehash;
                counter.rehash;
            }
        }
        cache[key] = value;
        counter[key] = 1;
    }
}

version(unittest) {
    /// keeps number of function evaluations
    static shared evaluations = 0;
    static shared int hits = 0;
    static shared int misses = 0;

    import std.stdio;
    import core.atomic;
    void printStats() {
        writeln("hits: ", hits, " misses: ", misses);
    }
}

auto memoize(alias func, uint maxItems=1024,
             alias CacheImpl=BasicCache, Args...)(Args args)
    if(isCallable!func && is(Args == ParameterTypeTuple!func)
                       && Args.length > 0) 
{
    alias ReturnType!func R; 

    static if (Args.length == 1) {
        static shared(CacheImpl!(maxItems, Args, R)) cache;
    } else {
        static shared(CacheImpl!(maxItems, Tuple!Args, R)) cache;
    }
    static shared bool init = false;

    if (!init) {
        synchronized {
            if (cache is null) {
                static if (Args.length == 1) {
                    cache = new shared(CacheImpl!(maxItems, Args, R))();
                } else {
                    cache = new shared(CacheImpl!(maxItems, Tuple!Args, R))();
                }
            }
            init = true;
        }
    }

    static if (Args.length == 1) {
        auto key = args;
    } else {
        auto key = tuple(args);
    }

    R* ret = (cast()cache).lookup(key);
    if (ret !is null) {
        version(unittest) {
            atomicOp!"+="(hits, 1);
        }
        return *ret;
    } else {
        version(unittest) {
            atomicOp!"+="(misses, 1);
        }
        synchronized(cache) {

            ret = (cast()cache).lookup(key);
            if (ret !is null) {
                return *ret;
            }

            version(unittest) {
                evaluations += 1;
            }

            auto result = func(args);
            (cast()cache).put(key, result);
            return result;
        }
    }
}

unittest {

    import core.thread;

    evaluations = 0;
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

    /// test FIFO cache
    alias memoize!(func, 2, FifoCache, int, int) fifomemoize;
    evaluations = 0;

    fifomemoize(1, 5); // 5
    fifomemoize(2, 3); // 5 6
    fifomemoize(2, 4); // 6 8
    fifomemoize(1, 5); // 8 5
    assert(evaluations == 4);
    fifomemoize(2, 4); // 8 5
    assert(evaluations == 4);
    fifomemoize(1, 7); // 5 7
    fifomemoize(1, 5); // 5 7
    assert(evaluations == 5);

    int foo(int x) {
        return x;
    }
    /// Test LU cache
    alias memoize!(foo, 3, LuCache, int) lumemoize;
    evaluations = 0;

    lumemoize(1);
    lumemoize(1);
    lumemoize(1); // 1
    lumemoize(2);
    lumemoize(2); // 1, 2
    lumemoize(3);
    lumemoize(3);
    lumemoize(3); // 1, 2, 3
    lumemoize(4); // 2 -> 4
    lumemoize(2); // 4 -> 2
    assert(evaluations == 5);
    lumemoize(3); // 1, 2, 3
    lumemoize(5); // 2 -> 5
    assert(evaluations == 6);
    lumemoize(4); // 5 -> 4
    lumemoize(9); // 4 -> 9
    assert(evaluations == 8);
}
