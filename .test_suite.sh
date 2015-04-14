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
    ./build/sambamba sort -t2 -m 50M ex1_header.nsorted.bam -o ex1_header.sorted.bam
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

testOverwriteProtection() {
    ./build/sambamba view ex1_header.nsorted.bam -f bam -o ./ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./build/sambamba merge ex1_header.sorted.bam ex1_header.sorted.bam ex1_header.sorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./build/sambamba sort -m 50M ex1_header.nsorted.bam -o ./build/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./build/sambamba markdup ex1_header.nsorted.bam ./build/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./build/sambamba slice ex1_header.nsorted.bam chr1 -o ./build/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
}

testSortingEmptyFile() {
    ./build/sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    ./build/sambamba sort -m 50M empty.bam -o empty2.bam 2>/dev/null
    assertEquals "0" `./build/sambamba view -c empty2.bam`
}

testMarkdupEmptyFile() {
    ./build/sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    ./build/sambamba markdup empty.bam empty.dedup.bam 2>/dev/null
    assertEquals "0" `./build/sambamba view -c empty.dedup.bam`
}

testCramView() {
    ./build/sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram
    ./build/sambamba view -C c1_pad2.cram >/dev/null
    assertEquals 0 $?
}

. shunit2-2.0.3/src/shell/shunit2
