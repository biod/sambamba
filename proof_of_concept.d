#!/usr/bin/rdmd
import std.parallelism;
import core.thread;
import std.stdio;
import std.datetime;

import rangetransformer;

/* imagine a situation: 

    you have to process the whole range using not much memory
    reading one block takes ~100ms
    while processing it takes ~300ms


That's the best we can hope for in terms of time:

time  input   thread   thread    thread
(ms)  range     1         2         3

  0    0        
100    1        0'
200    2        0'        1' 
300    3        0'-ready  1'        2'
400    4        3'        1'-ready  2'
500    5        3'        4'        2'-ready 
600    6        3'-ready  4'        5'
700    7        6'        4'-ready  5'
800    8        6'        7'        5'-ready
900    9        6'-ready  7'        8'
                9'        7'-ready  8'
                9'                  8'-ready
                9'-ready

New tasks are not run until old ones are removed via popFront(),
thus no more than 3 items are in task pool queue at every moment.

*/

enum N = 10;

class BgzfRange {

    auto delta() {
        return (Clock.currTime() - time).total!"msecs";
    }

    private uint i;

    SysTime time;

    this() {
        time = Clock.currTime();
        i = 0;
    }

    bool empty() @property {
        return i == N;
    }

    uint front() @property {
        return i;
    }

    void popFront() {
        writefln("getting element #%s | %sms", i, delta());
        Thread.sleep(dur!"msecs"(100));
        ++i;
    }
}

uint decompress(uint i) {
    writefln("Begin decompressing block %s", i);
    Thread.sleep(dur!"msecs"(300));
    writefln("End decompressing block %s", i);
    return i * i;
}

void main() {
    auto pool = new TaskPool(3);
    scope(exit) pool.finish();

    foreach (n; new RangeTransformer!(decompress, BgzfRange)(new BgzfRange(), pool)) {
        writeln(n);
    }
}
