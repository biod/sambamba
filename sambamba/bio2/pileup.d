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
module bio2.pileup;

import std.exception;
import std.stdio;

import std.experimental.logger;

import sambamba.bio2.constants;

immutable ulong DEFAULT_BUFFER_SIZE = 100_000;


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

alias ulong RingBufferIndex;

struct RingBuffer(T) {

  T[] _items;
  RingBufferIndex _head = 0;
  RingBufferIndex _tail = 0;

  /** initializes round buffer of size $(D n) */
  this(RingBufferIndex n) {
    _items = new T[n];
  }

  /*
  Does not work because data is no longer available!
  ~this() {
    // assert(is_empty); // make sure all items have been popped
  }
  */

  bool is_empty() @property const {
    return _tail == _head;
  }

  auto ref front() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[_head % $];
  }

  auto ref read_at(RingBufferIndex idx) {
    enforce(!is_empty, "ringbuffer is empty");
    enforce(idx >= _head, "ringbuffer range error");
    enforce(idx < _tail, "ringbuffer range error");
    return _items[idx % $];
  }

  void popFront() {
    enforce(!is_empty, "ringbuffer is empty");
    ++_head;
  }

  auto ref back() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[(_tail - 1) % $];
  }

  /// Puts item on the stack and returns the index
  RingBufferIndex put(T item) {
    enforce(!is_full, "ringbuffer is full - you need to expand buffer");
    enforce(_tail < _tail.max, "ringbuffer overflow");
    _items[_tail % $] = item; // uses copy semantics
    ++_tail;
    return _tail-1;
  }

  ulong length() @property const {
    writeln(_tail,":",_head,"= len ",_tail-_head);
    return _tail - _head;
  }

  bool is_full() @property const {
    return _items.length == length();
  }

  RingBufferIndex pushed() @property const {
    return _tail;
  }
  RingBufferIndex popped() @property const {
    return _head;
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

class PileUp(R) {
  RingBuffer!R ring;

  this(ulong bufsize=DEFAULT_BUFFER_SIZE) {
    ring = RingBuffer!R(bufsize);
  }

  RingBufferIndex push(R r) {
    return ring.put(r);
  }

  ref R front() {
    return ring.front();
  }

  ref R read_at(RingBufferIndex idx) {
    return ring.read_at(idx);
  }

  void popFront() {
    ring.popFront();
  }

  ulong depth(GenomePos pos, RingBufferIndex start_idx, RingBufferIndex stop_idx=RingBufferIndex.max) {
    size_t depth = 0;
    auto idx = start_idx;
    while (idx < ring._tail && idx < stop_idx) {
      auto read = ring._items[idx];
      writeln([idx,pos,read.start_pos,read.end_pos]);
      if (pos >= read.start_pos && pos <= read.end_pos) depth += 1;
      idx++;
    }
    return depth;
  }
}
