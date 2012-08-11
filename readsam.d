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
import samfile;
import std.stdio;
import std.parallelism;

void main(string[] args) {

    auto pool = new TaskPool(totalCPUs);
    scope(exit) pool.finish();
    auto sam = SamFile(args[1], pool);
    int i;
    foreach (read; sam.alignments) {
        if (read.read_name != "") {
            i += 1;
        }
    }
    writeln(i);
}
