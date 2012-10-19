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
__gshared ReferenceSequenceInfo[] refs;

DiploidCall5[] getSnps(P)(P pileup) {
    string r= "";
    if (pileup.ref_id != -1) {
        r = refs[pileup.ref_id].name;
    }
    return array((cast()caller).findSNPs(pileup, r));
}

void main(string[] args) {

    auto task_pool = new TaskPool(totalCPUs);
    scope(exit) task_pool.finish();

    auto fn = args[1];
    auto bam = BamFile(fn, task_pool);
    refs = bam.reference_sequences;
    auto reads = bam.alignments;

    caller = new shared(MaqSnpCaller)();
   
    // quality more than 50 is considered high
    (cast()caller).minimum_call_quality = 10.0f; 

    auto pileups = pileupChunks(reads, true, 4_000_000);

    foreach (snp; joiner(task_pool.map!getSnps(pileups, 32))) {
        if (snp.genotype.is_homozygous) {
            writeln(snp.position, '\t', snp.reference_base, '\t',
                    snp.genotype.base1, '\t', snp.quality);
        } else {
            writeln(snp.position, '\t', snp.reference_base, '\t',
                    snp.genotype.base1, ',', snp.genotype.base2, '\t', snp.quality);
        }
    }
}
