#!/usr/bin/env bash

testSortByName() {
    # use very tiny buffer of 200K so that multithreading is used
    ./build/sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam -m 200K
    ./build/sambamba view -t2 ex1_header.nsorted.bam > ex1_header.nsorted.sam
    assertEquals "3270" `wc -l < ex1_header.nsorted.sam`
    cat ex1_header.nsorted.sam | cut -f1 | LC_COLLATE=C sort -c
    assertEquals 0 $?
}

testSortByCoordinate() {
    ./build/sambamba sort -t2 ex1_header.nsorted.bam -o ex1_header.sorted.bam
    ./build/sambamba index -t2 ex1_header.sorted.bam
    assertEquals 0 $?
}

testSlice() {
    ./build/sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
    assertEquals 0 $?
}

testView() {
    assertEquals "1806" `./build/sambamba view -c ex1_header.sorted.bam chr2`
    assertEquals "1464" `./build/sambamba view -c ex1_header.sorted.bam chr1`
    assertEquals "0" `./build/sambamba view -c ex1_header.sorted.bam '*'`
}

. shunit2-2.0.3/src/shell/shunit2
