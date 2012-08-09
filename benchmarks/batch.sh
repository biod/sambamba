#!/bin/bash

# use some other file, no http/ftp support yet ;-)
BAMFILE=../../HG00476.chrom11.ILLUMINA.bwa.CHS.low_coverage.20111114.bam

FILENAME=benchmark.dat

echo -e "compiler\tn\tcompressed\tdecompressed\talignments" > $FILENAME

echo -ne "dmd\t" >> $FILENAME
./run_benchmarks $BAMFILE >> $FILENAME
echo -ne "gdc\t" >> $FILENAME
./run_benchmarks_gdc $BAMFILE >> $FILENAME

for i in {2..8}
do
    echo -ne "dmd\t" >> $FILENAME
    ./run_benchmarks_mt $BAMFILE $i >> $FILENAME
    echo -ne "gdc\t" >> $FILENAME
    ./run_benchmarks_mt_gdc $BAMFILE $i >> $FILENAME
done
