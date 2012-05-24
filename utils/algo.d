module utils.algo;

import std.range;

bool allDistinct(Range)(Range r) {
    uint items = 0;
    int[ElementType!Range] hash;
    foreach (elem; r) {
        items += 1;
        hash[elem] = 1;
    }
    return items == hash.keys.length;
}
