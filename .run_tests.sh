#!/usr/bin/env bash

testSortByName() {
    # use very tiny buffer of 200K so that multithreading is used
    ./build/sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam -m 200K
    ./build/sambamba view -t2 ex1_header.nsorted.bam > ex1_header.nsorted.sam
    assertEquals `wc -l < ex1_header.nsorted.sam` "3270"
    cat ex1_header.nsorted.sam | cut -f1 | LC_COLLATE=C sort -c
    assertEquals $? 0
}

testSortByCoordinate() {
    ./build/sambamba sort -t2 ex1_header.nsorted.bam -o ex1_header.sorted.bam
    ./build/sambamba index -t2 ex1_header.sorted.bam
    assertEquals $? 0
}

testSlice() {
    ./build/sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
    assertEquals $? 0
}

testView() {
    assertEquals `./build/sambamba view -c ex1_header.sorted.bam chr2` "1806"
    assertEquals `./build/sambamba view -c ex1_header.sorted.bam chr1` "1464"
    assertEquals `./build/sambamba view -c ex1_header.sorted.bam '*'` "0"
}

curl -L "http://downloads.sourceforge.net/shunit2/shunit2-2.0.3.tgz" | tar zx --overwrite
. shunit2-2.0.3/src/shell/shunit2
