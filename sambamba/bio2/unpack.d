/*
    This file is part of Sambamba.

    Copyright © 2012-2017    Pjotr Prins <pjotr.prins@thebird.nl>

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

/**
    Module for streaming unpacked (BAM) blocks
*/
module sambamba.bio2.unpack;

import std.stdio;

import sambamba.bio2.bgzf;

int unpack_bams(string[] infns, File outfile) {

  foreach (string fn; infns) {
    stderr.writeln(fn);

    foreach (immutable(ubyte[]) read; BgzfBlocks(fn)) {
      stdout.rawWrite(read);
    }
  }

  return 0;

}
