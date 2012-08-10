#!/bin/bash

BAMFILE=../../../HG00125.chrom20.ILLUMINA.bwa.GBR.low_coverage.20111114.bam

FILENAME=buffersize.dat

echo -e "size\tuser\tsystem\telapsed" > $FILENAME

for i in 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152
do
    echo -ne "$i\t" >> $FILENAME
    /usr/bin/time -f "%U\t%S\t%e" -o $FILENAME -a ./count $BAMFILE $i 
    sudo drop_cache.sh
done
