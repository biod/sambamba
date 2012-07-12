module common.nwayunion;

// copy-pasted from std.algorithm. The only difference is that front is not a ref.
import std.container;
import std.algorithm;
import std.functional;
import std.range;
import std.traits;

struct NWayUnion(alias less, RangeOfRanges)
{
    private alias .ElementType!(.ElementType!RangeOfRanges) ElementType;
    private alias binaryFun!less comp;
    private RangeOfRanges _ror;
    static bool compFront(.ElementType!RangeOfRanges a,
            .ElementType!RangeOfRanges b)
    {
        return comp(b.front, a.front);
    }
    BinaryHeap!(RangeOfRanges, compFront) _heap;

    this(RangeOfRanges ror)
    {
        _ror = remove!("a.empty", SwapStrategy.unstable)(ror);
        _heap.acquire(_ror);
    }

    @property bool empty() { return _ror.empty; }

    @property ElementType front() // <-------- the only difference
    {
        return _heap.front.front;
    }

    void popFront()
    {
        _heap.removeFront();
        _ror.back.popFront();
        if (_ror.back.empty)
        {
            _ror.popBack();
            return;
        }
        _heap.conditionalInsert(_ror.back) || assert(false);
    }
}

NWayUnion!(less, RangeOfRanges) nWayUnion
(alias less = "a < b", RangeOfRanges)
(RangeOfRanges ror)
{
    return typeof(return)(ror);
}
