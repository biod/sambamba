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
   Represent a pileup of reads in a buffer. When the tail reaches the
   buffer size the remainer gets copied to the start. Replacing it
   with a ringbuffer will be next (FIXME).
*/

class PileUp(R) {
  R[] buffer;
  size_t head = 0;
  size_t tail = 0;

  this(ulong bufsize=DEFAULT_BUFFER_SIZE) {
    buffer = new R[bufsize];
  }

  void push(R r) {
    enforce(tail < DEFAULT_BUFFER_SIZE);
    buffer[tail] = r;
    tail += 1;
    writeln("---->",tail-head);
  }

}
