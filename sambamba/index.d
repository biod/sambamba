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
module sambamba.index;

import std.stdio;
import std.stream;
import std.range;
import std.parallelism;
import std.getopt;

import sambamba.utils.common.progressbar;

import bio.bam.bai.indexing;
import bio.bam.reader;

void printUsage() {
    stderr.writeln("Usage: sambamba-index [-p|--show-progress] [-t|--nthreads NTHREADS] <input.bam> [<output.bai>]");
    stderr.writeln();
    stderr.writeln("\tIf output filename is not provided, appends '.bai' suffix");
    stderr.writeln("\tto the name of BAM file");
    stderr.writeln();
    stderr.writeln("\tIf -p option is specified, progressbar is shown in STDERR.");
    stderr.writeln("\tNumber of threads to use can be specified via -t option.");
}

version(standalone) {
    int main(string[] args) {
        return index_main(args);
    }
}

int index_main(string[] args) {

    bool show_progress;
    uint n_threads = totalCPUs;

    getopt(args,
           std.getopt.config.caseSensitive,
           "show-progress|p", &show_progress,
           "nthreads|t",      &n_threads);

    try {
        string out_filename = null;
        switch (args.length) {
            case 3:
                out_filename = args[2];
                goto case;
            case 2:
                if (out_filename is null)
                    out_filename = args[1] ~ ".bai";

                // default taskPool uses only totalCPUs-1 threads,
                // but in case of indexing the most time is spent
                // on decompression, and it makes perfect sense
                // to use all available cores for that
                //
                // (this is not the case with the sambamba tool where
                // filtering can consume significant amount of time)
                auto task_pool = new TaskPool(n_threads);
                scope(exit) task_pool.finish();

                auto bam = new BamReader(args[1], task_pool);
                Stream stream = new BufferedFile(out_filename, FileMode.Out);
                scope(exit) stream.close();

                if (show_progress) {
                    auto bar = new shared(ProgressBar)();
                    createIndex(bam, stream, (lazy float p) { bar.update(p); });
                    bar.finish();
                } else {
                    createIndex(bam, stream);
                }
                break;
            default:
                printUsage();
                return 0;
        }
    } catch (Throwable e) {
        stderr.writeln("sambamba-index: ", e.msg);
        version(development) { throw e; }
        return 1;
    }
    return 0;
}
