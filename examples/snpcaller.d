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
import std.getopt;

shared(MaqSnpCaller) caller;
__gshared ReferenceSequenceInfo[] refs;

DiploidCall5[] getSnps(P)(P pileup) {
    string r= "";
    if (pileup.ref_id != -1) {
        r = refs[pileup.ref_id].name;
    }
    return array((cast()caller).findSNPs(pileup, r));
}

void printUsage() {
    stderr.writeln(
"usage: ./snpcaller <input.bam> [OPTIONS]

OPTIONS:
    --out=FILENAME
        Write output to a file.

    --chr=CHROMOSOME_NAME
        Chromosome (if not given, only the first one presented in BAM file is processed).

    --start=START
        Position on the reference (0-based) from which to start calling. Default is 0.

    --stop=STOP
        Position on the reference (0-based) before which to stop calling. Default is -1UL.

    --min_call_quality=MIN_CALL_QUALITY
        Output only reads with quality higher than MIN_CALL_QUALITY. Default is 10. 

    --min_base_quality=MIN_BASE_QUALITY
        Discard reads with base quality less than MIN_BASE_QUALITY from computations.
        Default is 10.

    --threads=NTHREADS
        Number of threads to use. Default is the number of available cores plus one.

    --chunk_size=CHUNK_SIZE
        File is divided is chunks for parallelism. CHUNK_SIZE is approximate size in bytes
        of a single unpacked chunk. Thus, granularity is controlled by this parameter.
        Default value is 4000000 (4MB). Minimum is set internally to 10000.
");
}

void main(string[] args) {

    string chr = null;
    string output_filename = null;
    uint beg = 0;
    uint end = uint.max;
    int mincallquality = 10;
    ubyte minbasequality = 10;
    uint nthreads = totalCPUs + 1;
    size_t chunksize = 4_000_000;

    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               "out", &output_filename,
               "chr", &chr,
               "start", &beg,
               "stop", &end,
               "min_base_quality", &minbasequality,
               "min_call_quality", &mincallquality,
               "threads", &nthreads,
               "chunk_size", &chunksize);

        if (nthreads == 0) {
            nthreads = 1;
        }

        if (chunksize < 10_000) {
            chunksize = 10_000;
        }

    } catch (Exception e) {
        stderr.writeln(e.msg);
        stderr.writeln();
        printUsage();
        return;
    }

    if (args.length != 2) {
        printUsage();
        return;
    }

    try {

        auto task_pool = new TaskPool(nthreads - 1);
        scope(exit) task_pool.finish();

        auto fn = args[1];
        auto bam = BamFile(fn, task_pool);
        refs = bam.reference_sequences;

        File file = stdout;
        if (output_filename !is null) {
            file = File(output_filename, "w+");
        }

        caller = new shared(MaqSnpCaller)();

        // quality more than 50 is considered high
        (cast()caller).minimum_call_quality = cast(float)mincallquality;
        (cast()caller).minimum_base_quality = minbasequality;

        if (chr is null) {
            chr = refs[bam.alignments.front.ref_id].name;
        }

        auto reads = bam[chr][beg .. end];

        auto pileups = pileupChunks(reads, true, chunksize);

        foreach (snp; joiner(task_pool.map!getSnps(pileups, 32))) {
            if (snp.genotype.is_homozygous) {
                file.writeln(snp.position, '\t', snp.reference_base, '\t',
                        snp.genotype.base1, '\t', snp.quality);
            } else {
                file.writeln(snp.position, '\t', snp.reference_base, '\t',
                        snp.genotype.base1, ',', snp.genotype.base2, '\t', snp.quality);
            }

            file.flush();
        }

    } catch (Throwable e) {
        stderr.writeln(e.msg);
        return;
    }
}
