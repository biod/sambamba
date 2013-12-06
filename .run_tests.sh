#!/usr/bin/env bash
rdmd --compiler=ldmd2 -unittest -IBioD BioD/test/unittests.d

./build/sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam
./build/sambamba sort -t2 ex1_header.nsorted.bam -o ex1_header.sorted.bam
./build/sambamba index -t2 ex1_header.sorted.bam
./build/sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
./build/sambamba view -c ex1_header.sorted.bam chr2 > /dev/null
./build/sambamba view -c ex1_header.sorted.bam '*' > /dev/null
rm ex1_header*bam
