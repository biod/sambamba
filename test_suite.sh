#!/usr/bin/env bash

sambamba=./bin/sambamba
outdir=output
mkdir -p $outdir

nsortedbam=$outdir/ex1_header.nsorted.bam
sortedbam=$outdir/ex1_header.sorted.bam

testSamtoBam() {
    outfn=$nsortedbam
    $sambamba view -S test/ex1_header.sam -f bam > $outfn
    assertEquals "a51687bb6f8c0d5e7e8eecf08d84e2eb" $(md5sum $outfn |cut -c 1-32)
    $sambamba sort $outfn -o $sortedbam
    assertEquals "a51687bb6f8c0d5e7e8eecf08d84e2eb" $(md5sum $outfn |cut -c 1-32)
    $sambamba index $nsortedbam
    assertEquals "37a40998dbe8c74345005a4c343563ef" `./bin/sambamba view -f unpack $nsortedbam |md5sum|cut -c 1-32`
    assertEquals "6a48f9edf04bf95d3ac65d7395224440" `./bin/sambamba view -f unpack $sortedbam |md5sum|cut -c 1-32`
}

testSubSample() {
    # bam file is part of BioD
    outfn=$outdir/subsample.bam
    $sambamba subsample $outdir/ex1_header.sorted.bam -o$outfn > $outdir/subsample.out
    assertEquals "c3dc5e6d5dff1eaa3ff2d8e22cba5221" $(md5sum $outfn |cut -c 1-32)
    assertEquals "6a48f9edf04bf95d3ac65d7395224440" `./bin/sambamba view -f unpack output/subsample.bam |md5sum|cut -c 1-32`
}

testSortByName() {
    return
    # use very tiny buffer of 200K so that multithreading is used
    # $sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam -m 200K
    $sambamba view -t2 ex1_header.nsorted.bam > ex1_header.nsorted.sam
    assertEquals "3270" `wc -l < ex1_header.nsorted.sam`
    cat ex1_header.nsorted.sam | cut -f1 | LC_ALL=C sort -c
    assertEquals 0 $?
}

testSortByCoordinate() {
    return
    $sambamba sort -t2 -m 50M ex1_header.nsorted.bam -o ex1_header.sorted.bam
    $sambamba index -t2 ex1_header.sorted.bam
    assertEquals 0 $?
}

testSlice() {
    return
    $sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegions() {
    return
    $sambamba slice test/issue_204.bam 2:166860000-166870000 2:166860000-166870000 |\
      $sambamba view -c /dev/stdin > /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegionsBed() {
    return
    assertEquals "156" `$sambamba slice -L test/chr2_chr3_test_region.bed test/issue_204.bam |\
      $sambamba view -c /dev/stdin`
}

testView() {
    return
    assertEquals "1806" `$sambamba view -c ex1_header.sorted.bam chr2`
    assertEquals "1464" `$sambamba view -c ex1_header.sorted.bam chr1`
    assertEquals "0" `$sambamba view -c ex1_header.sorted.bam '*'`
}

testOverwriteProtection() {
    return
    $sambamba view ex1_header.nsorted.bam -f bam -o ./ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    $sambamba merge ex1_header.sorted.bam ex1_header.sorted.bam ex1_header.sorted.bam 2>/dev/null
    assertNotSame 0 $?
    $sambamba sort -m 50M ex1_header.nsorted.bam -o ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    $sambamba markdup ex1_header.nsorted.bam ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
    $sambamba slice ex1_header.nsorted.bam chr1 -o ./bin/../ex1_header.nsorted.bam 2>/dev/null
    assertNotSame 0 $?
}

testSortingEmptyFile() {
    return
    $sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    $sambamba sort -m 50M empty.bam -o empty2.bam 2>/dev/null
    assertEquals "0" `$sambamba view -c empty2.bam`
}

testMarkdupEmptyFile() {
    return
    $sambamba view ex1_header.sorted.bam -f bam -F "ref_id > 3" -o empty.bam 2>/dev/null
    $sambamba markdup empty.bam empty.dedup.bam 2>/dev/null
    assertEquals "0" `$sambamba view -c empty.dedup.bam`
}

testCramWriting() {
    return
    $sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram
    assertEquals 0 $?
}

testCramReading() {
    return
    $sambamba view -C c1_pad2.cram >/dev/null
    assertEquals 0 $?
}

testIndexUsage() {
    return
    rm *.bai && rm c1_*
    $sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o c1_pad2.cram &&
    $sambamba index -C c1_pad2.cram && test -e c1_pad2.cram.crai; assertEquals 0 $?
    $sambamba index -C c1_pad2.cram c1_cram_index && test -e c1_cram_index.crai; assertEquals 0 $?
    $sambamba index ex1_header.sorted.bam && test -e ex1_header.sorted.bam.bai; assertEquals 0 $?
    $sambamba index ex1_header.sorted.bam ex1_header.sorted.bai && test -e ex1_header.sorted.bai; assertEquals 0 $?
}

testIssue193() {
    return
    rm -f output_193.txt
    $sambamba depth base test/issue_193.bam > output_193.txt 2>/dev/null
    diff -q output_193.txt test/issue_193_expected_output.txt
    assertEquals 0 $?
}

testIssue204() {
    return
    rm -f output_204.txt
    $sambamba index test/issue_204.bam
    $sambamba depth region test/issue_204.bam -L 2:166868600-166868813 -T 15 -T 20 -T 25 -m > output_204.txt 2>/dev/null
    diff -q output_204.txt test/issue_204_expected_output.txt
    assertEquals 0 $?
}

testIssue206() {
    return
    $sambamba markdup ex1_header.sorted.bam ex1_header.dedup.bam 2>/dev/null
    $sambamba view -H ex1_header.dedup.bam | grep '@PG' | grep -q 'sambamba'
    assertEquals 0 $?

    $sambamba view ex1_header.sorted.bam -f bam -o ex1_header.filtered.bam\
                     -F "supplementary or secondary_alignment"
    $sambamba view -H ex1_header.filtered.bam | grep '@PG' | grep -q 'secondary'
    assertEquals 0 $?
}

testIssue225() {
    return
    $sambamba index test/issue225.bam
    $sambamba depth base -c 1 test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba depth base -c 0 test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?

    # exercise BED codepath as well
    $sambamba depth base -c 1 -L chrM test/issue225.bam > depth_base_225_1.txt 2>/dev/null
    diff -q depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba depth base -c 0 -L chrM test/issue225.bam > depth_base_225_0.txt 2>/dev/null
    diff -q depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?
}

shunit2=`which shunit2`
if [ -x "$shunit2" ]; then
    . $shunit2
else
    . shunit2-2.0.3/src/shell/shunit2
fi
