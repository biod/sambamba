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
module sambamba.index;

import std.stdio;
import undead.stream;
import std.range;
import std.parallelism;
import std.getopt;
import cram.reader;

import sambamba.utils.common.progressbar;

import bio.bam.bai.indexing;
import bio.bam.reader;

void printUsage() {
    stderr.writeln("Usage: sambamba-index [OPTIONS] <input.bam|input.cram> [output_file]");
    stderr.writeln();
    stderr.writeln("\tCreates index for a BAM or CRAM file");
    stderr.writeln();
    stderr.writeln("Options: -t, --nthreads=NTHREADS");
    stderr.writeln("               number of threads to use for decompression");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("               show progress bar in STDERR");
    stderr.writeln("         -c, --check-bins");
    stderr.writeln("               check that bins are set correctly");
    stderr.writeln("         -C, --cram-input");
    stderr.writeln("               specify that input is in CRAM format");
}

version(standalone) {
    int main(string[] args) {
        return index_main(args);
    }
}

int index_main(string[] args) {

    bool show_progress;
    bool check_bins;
    uint n_threads = totalCPUs;
    bool is_cram;
    string out_filename = null;

    getopt(args,
           std.getopt.config.caseSensitive,
           "show-progress|p", &show_progress,
           "nthreads|t",      &n_threads,
           "check-bins|c",    &check_bins,
           "cram-input|C",    &is_cram);

    try {
        if (args.length < 2 || args.length > 3) {
            printUsage();
            return 0;
        }

        if (!is_cram) {
            if (args.length > 2)
                out_filename = args[2];
            else
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
            bam.assumeSequentialProcessing();
            Stream stream = new BufferedFile(out_filename, FileMode.Out);
            scope(exit) stream.close();

            if (show_progress) {
                auto bar = new shared(ProgressBar)();
                createIndex(bam, stream, check_bins, (lazy float p) { bar.update(p); });
                bar.finish();
            } else {
                createIndex(bam, stream, check_bins);
            }
        } else {
            if (show_progress)
                stderr.writeln("[info] progressbar is unavailable for CRAM input");
            defaultPoolThreads = 0; // decompression not needed for CRAM
            auto cram = new CramReader(args[1], taskPool);
            cram.createIndex(args[$-1]);
        }
    } catch (Throwable e) {
        stderr.writeln("sambamba-index: ", e.msg);
        version(development) { throw e; }
        return 1;
    }
    return 0;
}
