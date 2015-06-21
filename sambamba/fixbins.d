/*
    This file is part of Sambamba.
    Copyright (C) 2012-2015    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.fixbins;

import std.stdio, std.getopt, std.parallelism;
import sambamba.utils.common.progressbar;
import sambamba.utils.common.overwrite;
import bio.bam.reader, bio.bam.writer;
import bio.bam.bai.bin;

void printUsage() {
    stderr.writeln("Usage: sambamba-fixbins [options] <input.bam> <output.bam>");
    stderr.writeln();
    stderr.writeln("Options: -t, --nthreads=NTHREADS");
    stderr.writeln("                    maximum number of threads to use");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("                    show progressbar in STDERR");
    stderr.writeln("         -l, --compression-level=LEVEL");
    stderr.writeln("                    specify compression level (from 0 to 9)");
}

version(standalone) {
    int main(string[] args) {
        return fixbins_main(args);
    }
}

int fixbins_main(string[] args) {
    try {
        int compression_level = -1;
        size_t nthreads = totalCPUs - 1;
        bool show_progress;
        getopt(args,
               std.getopt.config.caseSensitive,
               "nthreads|t", &nthreads,
               "compression-level|l", &compression_level,
               "show-progress|p", &show_progress);

        if (args.length != 3) {
            printUsage();
            return 0;
        }

        protectFromOverwrite(args[1], args[2]);
        
        auto bam = new BamReader(args[1]);
        bam.assumeSequentialProcessing();
        auto w = new BamWriter(args[2], compression_level);

        w.writeSamHeader(bam.header);
        w.writeReferenceSequenceInfo(bam.reference_sequences);

        void fixBins(R)(R reads) {
            foreach (r; reads) {
                auto start = r.position;
                auto end = start + r.basesCovered();
                uint correct_bin = cast(uint)reg2bin(start, end);
                uint* x = cast(uint*)r.raw_data.ptr + 2;
                *x = (*x & 0xFFFF) | (correct_bin << 16);
                w.writeRecord(r);
            }
        }

        if (show_progress) {
            auto bar = new shared(ProgressBar)();
            auto reads = bam.readsWithProgress((lazy float p) { bar.update(p); });
            fixBins(reads);
            bar.finish();
        } else {
            fixBins(bam.reads);
        }

        scope(exit) w.finish();

    } catch (Throwable e) {
        stderr.writeln("sambamba-fixbins: ", e.msg);
        version(development) { throw e; }
        return 1;
    }

    return 0;
}

