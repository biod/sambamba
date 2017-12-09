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

immutable ulong DEFAULT_BUFFER_SIZE = 1000_000;


/**
   Cyclic buffer or ringbuffer based on Artem's original. Uses copy
   semantics to copy a read in a pre-allocated buffer. New items get
   added to the tail, and used items get popped from the head
   (FIFO). Basically it is empty when pointers align and full when it
   is one before last. Not that the pointers are of size_t which puts
   a theoretical limit on the number of items that can be pushed.
*/

import core.stdc.string : memcpy;

struct RingBuffer(T) {

  private {
    T[] _items;
    size_t _put;
    size_t _taken;
  }

  /** initializes round buffer of size $(D n) */
  this(size_t n) {
    _items = new T[n];
  }

  ~this() {
    assert(is_empty); // make sure all items have been popped
  }

  // Input range primitives
  bool is_empty() @property const {
    return _put == _taken;
  }

  auto ref front() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[_taken % $];
  }

  void popFront() {
    enforce(!is_empty, "ringbuffer is empty");
    ++_taken;
  }

  auto ref back() @property {
    enforce(!is_empty, "ringbuffer is empty");
    return _items[(_put - 1) % $];
  }

  void put(T item) {
    enforce(!is_full, "ringbuffer is full");
    enforce(_put < _put.max, "ringbuffer overflow");
    // _items[_put % $] = item; // uses copy semantics
    auto b_item = _items[_put % $];
    memcpy(b_item.ptr,item.ptr,item.sizeof);
    ++_put;
  }

  bool is_full() @property const {
    return _put == _taken + _items.length;
  }

  /// Current number of elements
  size_t length() @property const {
    return _put - _taken;
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

  void push(R r) {
    ring.put(r);
    writeln("ring.length ",ring.length);
  }

}
