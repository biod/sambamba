module utils.algo;

bool allDistinct(Range)(Range r) {
    uint items = 0;
    int[typeof(r.front)] hash;
    foreach (elem; r) {
        items += 1;
        hash[elem] = 1;
    }
    return items == hash.keys.length;
}
