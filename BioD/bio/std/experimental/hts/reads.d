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

module bio.std.experimental.hts.reads;

import bio.std.experimental.hts.constants;
import bio.core.utils.exception;

import std.stdio;

bool read_overlaps(R)(GenomeLocation loc, R r) {
  assert(r.is_mapped);
  return r.ref_id == loc.ref_id && loc.pos >= r.start_pos && loc.pos <= r.end_pos;
}

bool reads_overlap(R)(R r1, R r2) {
  asserte(r1.is_mapped);
  asserte(r2.is_mapped);
  asserte(r2.ref_id == r1.ref_id);
  // r1                      rrrrrrrrrrr
  // r2                ---------???????????
  if (r2.start_pos < r1.start_pos && r2.end_pos >= r1.start_pos) {
    return true;
  }
  // r1                      rrrrrrrrrrr
  // r2                          ----?????
  return r2.start_pos >= r1.start_pos && r2.start_pos <= r1.end_pos;
}
