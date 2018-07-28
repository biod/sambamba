/*
    This file is part of Sambamba.
    Copyright (C) 2012-2017    Artem Tarasov <lomereiter@gmail.com>
    Copyright (C) 2012-2018    Pjotr Prins <pjotr.prins@thebird.nl>

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

import std.experimental.logger;

import sambamba.depth;
import sambamba.index;
import sambamba.fixbins;
import sambamba.flagstat;
import sambamba.markdup;
import sambamba.markdup2;
import sambamba.merge;
import sambamba.pileup;
import sambamba.sort;
import sambamba.slice;
import sambamba.subsample;
import sambamba.validate;
import sambamba.view;

import sambamba.utils.common.ldc_gc_workaround;

import utils.strip_bcf_header;
import utils.lz4;
import utils.version_ : VERSION;
import utils.ldc_version_info_ : LDC_VERSION_STRING, DMD_VERSION_STRING, LLVM_VERSION_STRING, BOOTSTRAP_VERSION_STRING;

import std.stdio;


void printUsage() {
    stderr.writeln("
Usage: sambamba [command] [args...]

  Available commands:

    view        view contents and convert from one format
                to another (SAM/BAM/CRAM/JSON/UNPACK)
    index       build index (BAI)
    merge       merge files (BAM)
    sort        sort file (BAM)
    slice       slice file (BAM using BED)
    markdup     mark or remove duplicates (BAM)
    subsample   subsamble (BAM)
    flagstat    output statistics (BAM)
    depth       output statistics (BAM)
    validate    simple validator (BAM)

  Work in progress (WIP):

    markdup2    mark or remove duplicates v2 (BAM)

  No longer recommended:

    mpileup     parallel execution of samtools (BAM)

To get help on a particular command, call it without args.

For bug reports and feature requests see

       https://github.com/biod/
");
}

void printVersion() {
    stderr.writeln();
    stderr.writeln("sambamba " ~ VERSION ~ " by Artem Tarasov and Pjotr Prins (C) 2012-2018");
    stderr.writeln("    LDC " ~ LDC_VERSION_STRING ~ " / DMD " ~ DMD_VERSION_STRING ~
     " / LLVM" ~ LLVM_VERSION_STRING ~ " / bootstrap " ~ BOOTSTRAP_VERSION_STRING);
    stderr.writeln();
}

int main(string[] args) {
    globalLogLevel(LogLevel.info);
    printVersion();

    if (args.length == 1) {
        printUsage();
        return 1;
    }

    auto _args = args[0] ~ args[2 .. $];

    switch (args[1]) {
        case "view":      return view_main(_args);
        case "index":     return index_main(_args);
        case "merge":     return merge_main(_args);
        case "sort":      return sort_main(_args);
        case "flagstat":  return flagstat_main(_args);
        case "slice":     return slice_main(_args);
        case "markdup":   return sambamba.markdup.markdup_main(_args);
        case "markdup2":  return sambamba.markdup2.markdup_main(_args);
        case "subsample": return subsample_main(_args);
        case "depth":     return depth_main(_args);
        case "mpileup":   return pileup_main(_args);
        case "validate":  return validate_main(_args);

        // hidden commands
        case "fixbins":   return fixbins_main(_args);
        case "strip_bcf_header": return strip_bcf_header_main(_args);
        case "lz4compress": return lz4compress_main();
        case "--version": printVersion(); return 0;
        default:
            printUsage();
            return 1;
    }
}
