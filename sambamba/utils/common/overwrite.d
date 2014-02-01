/*
    This file is part of Sambamba.
    Copyright (C) 2014    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.utils.common.overwrite;

import std.exception;
import std.path;

class OverwriteException : Exception {
    this(string output_fn) {
	super("Specified output filename " ~ output_fn ~
	      " clashes with one of input file names. Exiting.");
    }
}

void protectFromOverwrite(string input_filename, string output_filename) {
    auto input_fn = buildNormalizedPath(input_filename);
    auto output_fn = buildNormalizedPath(output_filename);
    if (input_fn.filenameCmp(output_fn) == 0)
	throw new OverwriteException(output_fn);
}
