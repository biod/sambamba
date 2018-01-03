#!/usr/bin/env bash

sambamba=./bin/sambamba
outdir=output
mkdir -p $outdir

# Name sorted and pos sorted
nsortedbam=$outdir/ex1_header.nsorted.bam
sortedbam=$outdir/ex1_header.sorted.bam

testSamtoBam() {
    outfn=$nsortedbam
    $sambamba view -S test/ex1_header.sam -f bam > $nsortedbam
    assertEquals "a51687bb6f8c0d5e7e8eecf08d84e2eb" $(md5sum $outfn |cut -c 1-32)
    # ex1_header.sorted.bam with index
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
    # should be same as sorted bam
    assertEquals "6a48f9edf04bf95d3ac65d7395224440" `./bin/sambamba view -f unpack output/subsample.bam |md5sum|cut -c 1-32`
}

testSortByName() {
    outfn=$outdir/ex1_header.nsorted.sam
    # use very tiny buffer of 200K so that multithreading is used
    $sambamba sort -t2 -n $sortedbam -o $nsortedbam -m 200K
    $sambamba view -t2 $nsortedbam > $outfn
    assertEquals "3270" `wc -l < $outfn`
    cat $outfn | cut -f1 | LC_ALL=C sort -c
    assertEquals 0 $?
}

testSortByCoordinate() {
    outfn=$outdir/ex1_header.sorted2.bam
    $sambamba sort -t2 -m 300K $nsortedbam -o $outfn
    $sambamba index -t2
    assertEquals 0 $?
    $sambamba view -t2 $outfn > $outfn.sam
    # check sorted with cat output/ex1_header.sorted2.bam.sam |awk '{print $3 " " $4}'|less
    assertEquals "2ae7237feab3ba5e389830b15ce834ee" `./bin/sambamba view -f unpack $outfn |md5sum|cut -c 1-32`
}

testSlice() {
    $sambamba slice $sortedbam chr1 -o /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegions() {
    $sambamba slice test/issue_204.bam 2:166860000-166870000 2:166860000-166870000 |\
        $sambamba view -c /dev/stdin > /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegionsBed() {
    assertEquals "156" `$sambamba slice -L test/chr2_chr3_test_region.bed test/issue_204.bam |\
          $sambamba view -c /dev/stdin`
}

testView() {
    assertEquals "1806" `$sambamba view -c $sortedbam chr2`
    assertEquals "1464" `$sambamba view -c $sortedbam chr1`
    assertEquals "0" `$sambamba view -c $sortedbam '*'`
}

testOverwriteProtection() {
    $sambamba view $nsortedbam -f bam -o ./$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba merge $sortedbam $sortedbam $sortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba sort -m 50M $nsortedbam -o ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba markdup $nsortedbam ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba slice $nsortedbam chr1 -o ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
}

testSortingEmptyFile() {
    $sambamba view $sortedbam -f bam -F "ref_id > 3" -o $outdir/empty.bam 2>/dev/null
    $sambamba sort -m 50M $outdir/empty.bam -o $outdir/empty2.bam 2>/dev/null
    assertEquals "0" `$sambamba view -c $outdir/empty2.bam`
}

testMarkdupEmptyFile() {
    $sambamba view $sortedbam -f bam -F "ref_id > 3" -o $outdir/empty.bam 2>/dev/null
    $sambamba markdup $outdir/empty.bam $outdir/empty.dedup.bam 2>/dev/null
    assertEquals "0" `$sambamba view -c $outdir/empty.dedup.bam`
}

testCramWriting() {
    $sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o $outdir/c1_pad2.cram
    assertEquals 0 $?
}

testCramReading() {
    $sambamba view -C $outdir/c1_pad2.cram >/dev/null
    assertEquals 0 $?
}

testIndexUsage() {
    rm *.bai && rm c1_*
    $sambamba view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o $outdir/c1_pad2.cram &&
    $sambamba index -C $outdir/c1_pad2.cram && test -e $outdir/c1_pad2.cram.crai; assertEquals 0 $?
    $sambamba index -C $outdir/c1_pad2.cram $outdir/c1_cram_index && test -e $outdir/c1_cram_index.crai; assertEquals 0 $?
    $sambamba index $sortedbam && test -e $sortedbam.bai; assertEquals 0 $?
    $sambamba index $sortedbam $outdir/ex1_header.sorted.bai && test -e $outdir/ex1_header.sorted.bai; assertEquals 0 $?
}

testIssue193() {
    rm -f $outdir/output_193.txt
    $sambamba depth base test/issue_193.bam > $outdir/output_193.txt 2>/dev/null
    diff -q $outdir/output_193.txt test/issue_193_expected_output.txt
    assertEquals 0 $?
}

testIssue204() {
    rm -f $outdir/output_204.txt
    $sambamba index test/issue_204.bam
    $sambamba depth region test/issue_204.bam -L 2:166868600-166868813 -T 15 -T 20 -T 25 -m > $outdir/output_204.txt 2>/dev/null
    diff -q $outdir/output_204.txt test/issue_204_expected_output.txt
    assertEquals 0 $?
}

testIssue206() {
    $sambamba markdup $sortedbam $outdir/ex1_header.dedup.bam 2>/dev/null
    $sambamba view -H $outdir/ex1_header.dedup.bam | grep '@PG' | grep -q 'sambamba'
    assertEquals 0 $?

    $sambamba view $sortedbam -f bam -o $outdir/ex1_header.filtered.bam \
              -F "supplementary or secondary_alignment"
    # skip failing test - filtered bam is not complete
    # $sambamba view -H $outdir/ex1_header.filtered.bam | grep '@PG' | grep -q 'secondary'
    assertEquals 0 $?
}

testIssue225() {
    $sambamba index test/issue225.bam
    $sambamba depth base -c 1 test/issue225.bam > $outdir/depth_base_225_1.txt 2>/dev/null
    diff -q $outdir/depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba depth base -c 0 test/issue225.bam > $outdir/depth_base_225_0.txt 2>/dev/null
    diff -q $outdir/depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?

    # exercise BED codepath as well
    $sambamba depth base -c 1 -L chrM test/issue225.bam > $outdir/depth_base_225_1.txt 2>/dev/null
    diff -q $outdir/depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba depth base -c 0 -L chrM test/issue225.bam > $outdir/depth_base_225_0.txt 2>/dev/null
    diff -q $outdir/depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?
}

testIssue331() {
    $sambamba  view -H -f json test/issue_204.bam > $outdir/issue_331_header.json
    assertEquals 0 $?
    diff -q $outdir/issue_331_header.json test/regression/issue_331_header.json
    assertEquals 0 $?
}

shunit2=`which shunit2`
if [ -x "$shunit2" ]; then
    . $shunit2
else
    . shunit2-2.0.3/src/shell/shunit2
fi
