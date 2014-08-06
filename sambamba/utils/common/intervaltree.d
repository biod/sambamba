// based on https://github.com/ekg/intervaltree
module sambamba.utils.common.intervaltree;

import std.range;
import std.algorithm;
import std.array;
import std.conv;

class IntervalTreeNode(T, pos_t) {
    pos_t start;
    pos_t stop;
    T value;

    this(pos_t start, pos_t stop, T value) {
        this.start = start;
        this.stop = stop;
        this.value = value;
    }

    override string toString() const {
        return "(" ~ start.to!string ~ " .. " ~ stop.to!string ~
               ", " ~ value.to!string ~ ")";
    }
}

class IntervalTree(T, pos_t) {
    alias IntervalTreeNode!(T, pos_t) interval;
    alias IntervalTree!(T, pos_t) intervalTree;

    interval[] intervals;
    intervalTree left;
    intervalTree right;
    pos_t center;

    this() {
        left = right = null;
        center = 0;
    }

    this(interval[] ivals, uint depth=16, uint minBucket=64,
         int leftExtent=0, int rightExtent=0)
    {
        left = right = null;
        intervals.length = ivals.length;
        intervals[] = ivals;

        --depth;
        if (depth == 0 || ivals.length < minBucket)
        {
            sort!((a, b) => a.start < b.start)(intervals[]);
        } else {
            if (leftExtent == 0 && rightExtent == 0) {
                sort!((a, b) => a.start < b.start)(intervals[]);
            }

            int leftP = 0;
            int rightP = 0;
            int centerP = 0;

            if (leftExtent != 0 || rightExtent != 0) {
                leftP = leftExtent;
                rightP = rightExtent;
            } else {
                leftP = intervals.front.start;
                rightP = intervals.map!(iv => iv.stop).reduce!max;
            }

            centerP = intervals[intervals.length / 2].start;
            center = centerP;

            interval[] lefts;
            interval[] rights;

            foreach (iv; ivals) {
                if (iv.stop <= center)
                    lefts ~= iv;
                else if (iv.start > center)
                    rights ~= iv;
                else intervals ~= iv;
            }

            if (lefts.length > 0)
                left = new intervalTree(lefts, depth, minBucket,
                                        leftP, centerP);

            if (rights.length > 0)
                right = new intervalTree(rights, depth, minBucket,
                                         centerP, rightP);
        }
    }

    void findOverlapping(pos_t start, pos_t stop,
                               scope void delegate(ref interval) cb)
    {
        if (!intervals.empty && !(stop <= intervals.front.start)) {
            foreach (ref iv; intervals)
                if (iv.stop > start && iv.start < stop)
                    cb(iv);

            if (left !is null && start < center)
                left.findOverlapping(start, stop, cb);

            if (right !is null && stop >= center)
                right.findOverlapping(start, stop, cb);
        }
    }
}
