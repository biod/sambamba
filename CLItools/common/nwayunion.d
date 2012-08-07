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
