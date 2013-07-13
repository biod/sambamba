/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
import sambamba.view;
import sambamba.index;
import sambamba.merge;
import sambamba.sort;
import sambamba.flagstat;
import sambamba.slice;

import sambamba.utils.common.ldc_gc_workaround;

import std.stdio;

void printUsage() {
    stderr.writeln("sambamba v0.3.3");
    stderr.writeln();
    stderr.writeln("Usage: sambamba [command] [args...]");
    stderr.writeln();
    stderr.writeln("    Available commands: 'view', 'index', 'merge', 'sort', 'flagstat', 'slice'");
    stderr.writeln("    To get help on a particular command, just call it without args.");
}

int main(string[] args) {
    if (args.length == 1) {
        printUsage();
        return 1;
    }

    auto _args = args[0] ~ args[2 .. $];

    switch (args[1]) {
        case "view":     return view_main(_args);
        case "index":    return index_main(_args);
        case "merge":    return merge_main(_args);
        case "sort":     return sort_main(_args);
        case "flagstat": return flagstat_main(_args);
        case "slice":    return slice_main(_args);
        default: 
            printUsage();
            return 1;
    }
}
