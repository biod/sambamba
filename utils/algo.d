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

import std.functional;

auto argmax(alias func, S)(S set) {
    auto best_elem = set.front;
    auto best_value = unaryFun!func(best_elem);

    set.popFront();
    foreach (elem; set) {
        auto value = unaryFun!func(elem);
        if (value > best_value) {
            best_value = value;
            best_elem = elem;
        }
    }

    return best_elem;
}
