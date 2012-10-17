import bamfile;
import snpcallers.maq;
import pileuprange;

import splitter;

import BioD.Call;

import std.stdio;
import std.conv;
import std.array;
import std.algorithm;
import std.parallelism;
import std.range;

shared(MaqSnpCaller) caller;

DiploidCall5[] getSnps(P)(P pileup) {
    return array((cast()caller).findSNPs(pileup));
}

void main(string[] args) {

    auto task_pool = new TaskPool(totalCPUs);
    scope(exit) task_pool.finish();

    auto fn = args[1];
    auto bam = BamFile(fn, task_pool);
    auto reads = bam.alignments;

    caller = new shared(MaqSnpCaller)();

    auto pileups = pileupChunks(reads, 8_000_000);

    foreach (snp; joiner(task_pool.map!getSnps(pileups, 32))) {
        writeln(snp.position, " ", snp.reference_base, " ", snp.genotype, " ", snp.quality);
    }
}
