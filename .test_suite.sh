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

testIndexUsage() {
    rm *.bai && rm c1_*
    ./build/sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram &&
    ./build/sambamba index -C c1_pad2.cram && test -e c1_pad2.cram.crai
    ./build/sambamba index -C c1_pad2.cram c1_cram_index && test -e c1_cram_index.crai
    ./build/sambamba index ex1_header.sorted.bam && test -e ex1_header.sorted.bam.bai
    ./build/sambamba index ex1_header.sorted.bam ex1_header.sorted.bai && test -e ex1_header.sorted.bai
    assertEquals 0 $?
}

testIssue193() {
    rm -f output_193.txt
    ./build/sambamba depth base test/issue_193.bam > output_193.txt 2>/dev/null
    diff -q output_193.txt test/issue_193_expected_output.txt
    assertEquals 0 $?
}

testIssue204() {
    rm -f output_204.txt
    ./build/sambamba index test/issue_204.bam
    ./build/sambamba depth region test/issue_204.bam -L 2:166868600-166868813 -T 15 -T 20 -T 25 -m > output_204.txt 2>/dev/null
    diff -q output_204.txt test/issue_204_expected_output.txt
    assertEquals 0 $?
}

testIssue206() {
    ./build/sambamba markdup ex1_header.sorted.bam ex1_header.dedup.bam 2>/dev/null
    ./build/sambamba view -H ex1_header.dedup.bam | grep '@PG' | grep -q 'sambamba markdup'
    assertEquals 0 $?
}

testIssue225() {
    ./build/sambamba index test/issue225.bam
    ./build/sambamba depth base -c 1 test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    ./build/sambamba depth base -c 0 test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?

    # exercise BED codepath as well
    ./build/sambamba depth base -c 1 -L chrM test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    ./build/sambamba depth base -c 0 -L chrM test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?
}

. shunit2-2.0.3/src/shell/shunit2
