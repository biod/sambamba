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
import contrib.undead.stream;
import std.range;
import std.parallelism;
import std.getopt;
import std.file;

import sambamba.utils.common.progressbar;

import bio.std.hts.bam.bai.indexing;
import bio.std.hts.bam.reader;
import bio.std.file.fai;

void printUsage() {
    stderr.writeln("Usage: sambamba-index [OPTIONS] <input.bam|input.fasta> [output_file]");
    stderr.writeln();
    stderr.writeln("\tCreates index for a BAM, or FASTA file");
    stderr.writeln();
    stderr.writeln("Options: -t, --nthreads=NTHREADS");
    stderr.writeln("               number of threads to use for decompression");
    stderr.writeln("         -p, --show-progress");
    stderr.writeln("               show progress bar in STDERR");
    stderr.writeln("         -c, --check-bins");
    stderr.writeln("               check that bins are set correctly");
    stderr.writeln("         -F, --fasta-input");
    stderr.writeln("               specify that input is in FASTA format");
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
    bool is_fasta;
    string out_filename = null;

    // Getopt removes parsed arguments
    getopt(args,
           std.getopt.config.caseSensitive,
           "show-progress|p", &show_progress,
           "nthreads|t",      &n_threads,
           "check-bins|c",    &check_bins,
           "fasta-input|F",    &is_fasta);

    try {
        // args = ["sambamba", "input.bam", "output.bai"]
        if (args.length != 2 && args.length != 3) {
            printUsage();
            return 0;
        }

        string input_filename = args[1];

        if (!is_cram && !is_fasta) {
            if (args.length > 2)
                out_filename = args[2];
            else
                out_filename = input_filename ~ ".bai";

            // default taskPool uses only totalCPUs-1 threads,
            // but in case of indexing the most time is spent
            // on decompression, and it makes perfect sense
            // to use all available cores for that
            //
            // (this is not the case with the sambamba tool where
            // filtering can consume significant amount of time)
            auto task_pool = new TaskPool(n_threads);
            scope(exit) task_pool.finish();

            auto bam = new BamReader(input_filename, task_pool);
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
        }
        else if(is_fasta) {
            stderr.writeln("Indexing FASTA file...");
            if (show_progress)
                stderr.writeln("[info] progressbar is unavailable for FASTA input");
            if (args.length > 2)
                out_filename = args[2];
            else
                out_filename = input_filename ~ ".fai";

            string records;
            foreach(FaiRecord rec; buildFai(input_filename))
                records ~= rec.toString() ~ '\n';

            std.file.write(out_filename, records);

        }
    } catch (Throwable e) {
        stderr.writeln("sambamba-index: ", e.msg);
        version(development) { throw e; }
        return 1;
    }
    return 0;
}
