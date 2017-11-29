#!/usr/bin/env bash

testSortByName() {
    # use very tiny buffer of 200K so that multithreading is used
    # ./bin/sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam -m 200K
    ./bin/sambamba view -t2 ex1_header.nsorted.bam > ex1_header.nsorted.sam
    assertEquals "3270" `wc -l < ex1_header.nsorted.sam`
    cat ex1_header.nsorted.sam | cut -f1 | LC_ALL=C sort -c
    assertEquals 0 $?
}

testSortByCoordinate() {
    ./bin/sambamba sort -t2 -m 50M ex1_header.nsorted.bam -o ex1_header.sorted.bam
    ./bin/sambamba index -t2 ex1_header.sorted.bam
    assertEquals 0 $?
}

testSlice() {
    ./bin/sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegions() {
    ./bin/sambamba slice test/issue_204.bam 2:166860000-166870000 2:166860000-166870000 |\
      ./bin/sambamba view -c /dev/stdin > /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegionsBed() {
    assertEquals "156" `./bin/sambamba slice -L test/chr2_chr3_test_region.bed test/issue_204.bam |\
      ./bin/sambamba view -c /dev/stdin`
}

testView() {
    assertEquals "1806" `./bin/sambamba view -c ex1_header.sorted.bam chr2`
    assertEquals "1464" `./bin/sambamba view -c ex1_header.sorted.bam chr1`
    assertEquals "0" `./bin/sambamba view -c ex1_header.sorted.bam '*'`
}

testOverwriteProtection() {
    ./bin/sambamba view ex1_header.nsorted.bam -f bam -o ./ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./bin/sambamba merge ex1_header.sorted.bam ex1_header.sorted.bam ex1_header.sorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./bin/sambamba sort -m 50M ex1_header.nsorted.bam -o ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./bin/sambamba markdup ex1_header.nsorted.bam ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    ./bin/sambamba slice ex1_header.nsorted.bam chr1 -o ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
}

testSortingEmptyFile() {
    ./bin/sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    ./bin/sambamba sort -m 50M empty.bam -o empty2.bam 2>/dev/null
    assertEquals "0" `./bin/sambamba view -c empty2.bam`
}

testMarkdupEmptyFile() {
    ./bin/sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    ./bin/sambamba markdup empty.bam empty.dedup.bam 2>/dev/null
    assertEquals "0" `./bin/sambamba view -c empty.dedup.bam`
}

testCramWriting() {
    ./bin/sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram
    assertEquals 0 $?
}

testCramReading() {
    ./bin/sambamba view -C c1_pad2.cram >/dev/null
    assertEquals 0 $?
}

testIndexUsage() {
    rm *.bai && rm c1_*
    ./bin/sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram &&
    ./bin/sambamba index -C c1_pad2.cram && test -e c1_pad2.cram.crai; assertEquals 0 $?
    ./bin/sambamba index -C c1_pad2.cram c1_cram_index && test -e c1_cram_index.crai; assertEquals 0 $?
    ./bin/sambamba index ex1_header.sorted.bam && test -e ex1_header.sorted.bam.bai; assertEquals 0 $?
    ./bin/sambamba index ex1_header.sorted.bam ex1_header.sorted.bai && test -e ex1_header.sorted.bai; assertEquals 0 $?
}

testIssue193() {
    rm -f output_193.txt
    ./bin/sambamba depth base test/issue_193.bam > output_193.txt 2>/dev/null
    diff -q output_193.txt test/issue_193_expected_output.txt
    assertEquals 0 $?
}

testIssue204() {
    rm -f output_204.txt
    ./bin/sambamba index test/issue_204.bam
    ./bin/sambamba depth region test/issue_204.bam -L 2:166868600-166868813 -T 15 -T 20 -T 25 -m > output_204.txt 2>/dev/null
    diff -q output_204.txt test/issue_204_expected_output.txt
    assertEquals 0 $?
}

testIssue206() {
    ./bin/sambamba markdup ex1_header.sorted.bam ex1_header.dedup.bam 2>/dev/null
    ./bin/sambamba view -H ex1_header.dedup.bam | grep '@PG' | grep -q 'sambamba'
    assertEquals 0 $?

    ./bin/sambamba view ex1_header.sorted.bam -f bam -o ex1_header.filtered.bam\
                     -F "supplementary or secondary_alignment"
    ./bin/sambamba view -H ex1_header.filtered.bam | grep '@PG' | grep -q 'secondary'
    assertEquals 0 $?
}

testIssue225() {
    ./bin/sambamba index test/issue225.bam
    ./bin/sambamba depth base -c 1 test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    ./bin/sambamba depth base -c 0 test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?

    # exercise BED codepath as well
    ./bin/sambamba depth base -c 1 -L chrM test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    ./bin/sambamba depth base -c 0 -L chrM test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?
}

shunit2=`which shunit2`
if [ -x "$shunit2" ]; then
    . $shunit2
else
    . shunit2-2.0.3/src/shell/shunit2
fi
