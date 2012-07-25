module sambamba.index;

import std.stdio;
import std.stream;
import std.range;
import std.parallelism;
import std.getopt;

import common.progressbar;

import bai.indexing;
import bamfile;

void printUsage() {
    writeln("Usage: sambamba-index [-p|--show-progress] <input.bam> [<output.bai>]");
    writeln();
    writeln("\tIf output filename is not provided, appends '.bai' suffix");
    writeln("\tto the name of BAM file");
    writeln();
    writeln("\tIf -p option is specified, progressbar is shown in STDERR.");
}

version(standalone) {
    int main(string[] args) {
        return index_main(args);
    }
}

int index_main(string[] args) {

    bool show_progress;

    getopt(args,
           std.getopt.config.caseSensitive,
           "show-progress|p", &show_progress);

    try {
        string out_filename = null;
        switch (args.length) {
            case 3:
                out_filename = args[2];
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
                auto task_pool = new TaskPool(totalCPUs);
                scope(exit) task_pool.finish();

                auto bam = BamFile(args[1], task_pool);
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
        return 1;
    }
    return 0;
}
