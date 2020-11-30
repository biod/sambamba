/*
    This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/
module bio.std.experimental.hts.pileup;

import std.conv;
import std.exception;
import std.stdio;
import std.traits;
import std.typecons;

import std.experimental.logger;

import bio.std.experimental.hts.constants;
import bio.core.utils.exception;

immutable ulong DEFAULT_BUFFER_SIZE = 1_000_000;

/**
   Cyclic buffer or ringbuffer based on Artem's original. Uses copy
   semantics to copy a read in a pre-allocated buffer. New items get
   added to the tail, and used items get popped from the head
   (FIFO). Basically it is empty when pointers align and full when
   head - tail equals length. Not that the pointers are of size_t
   which puts a theoretical limit on the number of items that can be
   pushed.
*/

import core.stdc.string : memcpy;

// alias ulong RingBufferIndex;

struct RingBufferIndex {
  alias Representation = ulong;
  private ulong value = 0;

  this(ulong v) {
    value = v;
  }

  // @disable this(this); // disable copy semantics;

  auto get() inout {
    return value;
  }

  auto max() @property {
    return value.max;
  }

  void opAssign(U)(U rhs) if (is(typeof(Checked!(T, Hook)(rhs)))) {
    value = rhs;
  }

  bool opEquals(U, this _)(U rhs) {
    return value == rhs;
  }

  auto opCmp(U, this _)(const U rhs) {
    return value < rhs ? -1 : value > rhs;
  }

  ulong opUnary(string s)() if (s == "++") {
    return ++value;
  }

}


struct RingBuffer(T) {

  T[] _items;
  RingBufferIndex _head;
  RingBufferIndex _tail;
  size_t max_size = 0;

  /** initializes round buffer of size $(D n) */
  this(size_t n) {
    _items = new T[n];
    // _items.reserve(n);
  }

  @disable this(this); // disable copy semantics;

  /*
  Does not work because data is no longer available!
  ~this() {
    // assert(is_empty); // make sure all items have been popped
  }
  */

  bool empty() @property @nogc nothrow const {
    return _tail == _head;
  }

  alias empty is_empty;

  auto ref front() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[_head.get() % $];
  }
  alias back last;

  auto ref back() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[(_tail.get() - 1) % $];
  }

  alias front first;

  bool is_tail(RingBufferIndex idx) {
    return idx == _tail.get()-1;
  }

  ref T get_at(RingBufferIndex idx) {
    enforce(!is_empty, "ringbuffer is empty");
    enforce(idx >= _head, "ringbuffer range error (idx before front)");
    enforce(idx != _tail, "ringbuffer range error (idx at end)");
    enforce(idx < _tail, "ringbuffer range error (idx after end)");
    return _items[idx.get() % $];
  }

  bool is_valid(RingBufferIndex idx) {
    enforce(!is_empty, "ringbuffer is empty");
    enforce(idx >= _head, "ringbuffer range error (idx before front)");
    enforce(idx != _tail, "ringbuffer range error (idx at end)");
    enforce(idx < _tail, "ringbuffer range error (idx after end)");
    return true;
  }

  // This function is a hack.
  void update_at(RingBufferIndex idx, T item) {
    is_valid(idx);
    _items[idx.get() % $] = item; // uses copy semantics
  }

  RingBufferIndex popFront() {
    enforce(!is_empty, "ringbuffer is empty");
    static if (__traits(compiles, _items[0].cleanupx)) {
      // write("x");
      _items[_head.get() % $].cleanup();
    }
    ++_head.value;
    return _head;
  }

  /// Puts item on the stack and returns the index
  RingBufferIndex put(T item) {
    enforce(!is_full, "ringbuffer is full - you need to expand buffer");
    enforce(_tail < _tail.max, "ringbuffer overflow");
    max_size = length > max_size ? length : max_size;
    _items[_tail.get() % $] = item; // uses copy semantics
    auto prev = _tail;
    ++_tail.value;
    return prev;
  }

  ulong length() @property const {
    // writeln(_tail.get(),":",_head.get(),"= len ",_tail.get()-_head.get());
    return _tail.get() - _head.get();
  }

  bool is_full() @property const {
    return _items.length == length();
  }

  bool in_range(RingBufferIndex idx) @property const {
    return idx >= _head && idx < _tail;
  }

  ulong pushed() @property const {
    return _tail.value;
  }
  ulong popped() @property const {
    return _head.value;
  }

  @property void cleanup() {
    _head = RingBufferIndex();
    _tail = RingBufferIndex();
  }

  string toString() {
    string res = "ring ";
    for(RingBufferIndex i = _head; i<_tail; i++)
      res ~= to!string(get_at(i));
    return res;
  }

  @property string stats() {
    return "Ringbuffer pushed " ~ to!string(pushed) ~ " popped " ~ to!string(popped) ~ " max-size " ~
      to!string(max_size) // , "/", (pileup.ring.max_size+1)/pileup.ring.length);
      ;
  }
}

unittest {
    auto buf = RingBuffer!int(4);
    assert(buf.is_empty);

    buf.put(1);
    buf.put(2);
    assert(buf.length == 2);
    assert(buf.front == 1);
    buf.popFront(); // 1
    buf.popFront(); // 2
    buf.put(2);
    buf.put(1);
    buf.put(0);
    buf.put(3);
    assert(buf.is_full);
    assert(buf.front == 2);
    buf.popFront();
    assert(buf.front == 1);
    buf.put(4);
    buf.popFront();
    assert(buf.front == 0);
    buf.popFront();
    assert(buf.front == 3);
    buf.popFront();
    assert(buf.front == 4);
    buf.popFront();
    assert(buf.is_empty);
}

/**
   Represent a pileup of reads in a buffer.
*/

/*
  Read RingBuffer with current pointer, so you have three states
  (first, current, last).
*/

class PileUp(R) {
  RingBuffer!R ring;
  Nullable!RingBufferIndex current;

  this(ulong bufsize=DEFAULT_BUFFER_SIZE) {
    ring = RingBuffer!R(bufsize);
    set_current_to_head;
  }

  RingBufferIndex push(R r) { return ring.put(r); }
  bool empty() @property const { return ring.is_empty();}
  bool is_full() @property const { return ring.is_full();}
  RingBufferIndex popFront() { return ring.popFront(); }
  ref R front() { return ring.front(); }
  alias front leftmost;
  ref R rightmost() { return ring.back(); }
  ref R read(RingBufferIndex idx) {
    enforce(ring.in_range(idx), "idx should be set for PileUp.read");
    return ring.get_at(idx);
  }
  ref R read_current() {
    enforce(!current.isNull, "current should be set for PileUp.read_current");
    return read(current);
  }
  bool is_at_end(RingBufferIndex idx) { return ring.is_tail(idx); }

  @property void current_inc() {
    asserte(!empty);
    asserte(!ring.is_tail(current));
    ++current;
  }

  @property void set_current_to_head() {
    current = ring._head; // note pileup can be empty
  }

  void current_reset() {
    current = RingBufferIndex();
  }

  @property bool current_is_tail() {
    return ring.is_tail(current);
  }

  void each(void delegate(R) dg) {
    auto idx = ring._head;
    while(!ring.is_tail(idx)) {
      auto r = read(idx);
      dg(r);
      idx++;
    }
  }

  void each_left_of_current(void delegate(RingBufferIndex, R) dg) {
    R cur = read_current;
    if (cur.is_unmapped) return;
    auto idx = ring._head;
    while(!ring.is_tail(idx)) {
      auto r = read(idx);
      if (r.is_mapped && r.end_pos >= cur.start_pos)
        return;
      dg(idx,r);
      idx++;
    }
  }

  void purge_while(bool delegate(R) dg) {
    while(!empty) {
      if (!dg(front))
        return; // skip the rest
      popFront();
    }
  }

  void purge(void delegate(R) dg) {
    while(!empty) {
      dg(front);
      popFront();
    }
    set_current_to_head();
  }

  @property string stats() {
    return ring.stats();
  }

  override string toString() {
    return stats;
  }
}
