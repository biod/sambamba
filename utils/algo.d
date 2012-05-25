module utils.algo;

import std.range;

/**
  This function is supposed to be used on a small amount of objects,
  typically tags of the same alignment. Therefore it runs in O(N^2)
  but doesn't require any memory allocations.
*/
bool allDistinct(Range)(Range r) {
    uint sz = 0;
    uint eq = 0;
    foreach (e1; r) {
        ++sz;
        foreach (e2; r) {
            if (e1 == e2) {
                eq += 1;
            }
        }
    }
    return sz == eq;
}

private import std.algorithm;
static if (!__traits(compiles, any!"a == 2"([1,2,3]))) {
    /** GDC uses older phobos library, so let's define 'all' and 'any' functions */

    import std.functional;

    /// check if all elements satisfy the condition
    bool all(alias F, R)(R r) {
        foreach (e; r) {
            if (!unaryFun!F(e)) {
                return false;
            }
        }
        return true;
    }

    /// check if any element satisfies the condition
    bool any(alias F, R)(R r) {
        foreach (e; r) {
            if (unaryFun!F(e)) {
                return true;
            }
        }
        return false;
    }
}
