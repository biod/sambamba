// based on https://github.com/ekg/intervaltree
module sambamba.utils.common.intervaltree;

import std.range;
import std.algorithm;
import std.array;
import std.conv;

class IntervalTreeNode(T, pos_t) {
    static if (is(T == void)) {
        pos_t start;
        pos_t stop;

        this(pos_t start, pos_t stop) {
            this.start = start;
            this.stop = stop;
        }

        override string toString() const {
            return "(" ~ start.to!string ~ " .. " ~ stop.to!string ~ ")";
        }
    } else {
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
}

class IntervalTree(T, pos_t) {
    alias IntervalTreeNode!(T, pos_t) interval;
    alias IntervalTree!(T, pos_t) intervalTree;

    void print(int depth=0) {
        import std.stdio;
        writeln(repeat(' ', depth), "INTERVALS: ", intervals);
        if (left !is null) {
            writeln(repeat(' ', depth), "LEFT");
            left.print(depth + 4);
        }
        if (right !is null) {
            writeln(repeat(' ', depth), "RIGHT");
            right.print(depth + 4);
        }
    }

    interval[] intervals;
    intervalTree left;
    intervalTree right;
    pos_t center;

    this() {
        left = right = null;
        center = 0;
    }

    this(ref interval[] ivals, uint depth=16, uint minBucket=2,//64,
         int leftExtent=0, int rightExtent=0)
    {
        left = right = null;

        --depth;
        if (depth == 0 || ivals.length < minBucket)
        {
            sort!((a, b) => a.start < b.start)(ivals[]);
            intervals = ivals.dup;
        } else {
            if (leftExtent == 0 && rightExtent == 0) {
                sort!((a, b) => a.start < b.start)(ivals[]);
            }

            int leftP = 0;
            int rightP = 0;
            int centerP = 0;

            if (leftExtent != 0 || rightExtent != 0) {
                leftP = leftExtent;
                rightP = rightExtent;
            } else {
                leftP = ivals.front.start;
                rightP = ivals.map!(iv => iv.stop).reduce!max;
            }

            centerP = ivals[ivals.length / 2].start;
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

    static struct OverlapWalker {
        pos_t start;
        pos_t stop;
        intervalTree tree;
        int opApply(scope int delegate(ref interval) dg) {
            int result = 0;
            if (!tree.intervals.empty && !(stop <= tree.intervals.front.start)) {
                foreach (ref iv; tree.intervals)
                    if (iv.stop > start && iv.start < stop) {
                        result = dg(iv);
                        if (result != 0)
                            return result;
                    }
            }

            if (tree.left !is null && start < tree.center) {
                result = tree.left.eachOverlap(start, stop).opApply(dg);
                if (result != 0)
                    return result;
            }

            if (tree.right !is null && stop >= tree.center) {
                result = tree.right.eachOverlap(start, stop).opApply(dg);
                if (result != 0)
                    return result;
            }

            return result;
        }
    }

    auto eachOverlap(pos_t start=pos_t.min, pos_t stop=pos_t.max)
    {
        return OverlapWalker(start, stop, this);
    }
}
