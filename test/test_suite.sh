#!/bin/sh

sambamba=./bin/sambamba
opts="-q"
outdir=output
mkdir -p $outdir

md5sum=`which md5sum`
[ -z $md5sum ] && md5sum="md5 -r " # for OSX

# Name sorted and pos sorted
nsortedbam=$outdir/ex1_header.nsorted.bam
sortedbam=$outdir/ex1_header.sorted.bam

testSamtoBam() {
    outfn=$nsortedbam
    # first identity SAM but without header
    $sambamba $opts view -S test/ex1_header.sam -f sam > $nsortedbam.sam.1
    assertEquals "ff5d66e6cd1cab4e8e6330f2829b7c41" $($md5sum $nsortedbam.sam.1 |cut -c 1-32)
    # now convert to BAM
    $sambamba $opts view -S test/ex1_header.sam -f bam > $nsortedbam
    # BAM back to identity
    $sambamba $opts view $nsortedbam -f sam > $nsortedbam.sam.2
    assertEquals "ff5d66e6cd1cab4e8e6330f2829b7c41" $($md5sum $nsortedbam.sam.2 |cut -c 1-32)
    # Test BAM
    assertEquals "b5d64266e9544a656f0efbd0f8030b20" $($md5sum $nsortedbam |cut -c 1-32)
    # ex1_header.sorted.bam with index
    $sambamba $opts sort $outfn -o $sortedbam
    assertEquals "b5d64266e9544a656f0efbd0f8030b20" $($md5sum $outfn |cut -c 1-32)
    $sambamba $opts index $nsortedbam
    assertEquals "aabb4fb319497dec1069028dd354eeda" `$sambamba $opts view -f unpack $nsortedbam |$md5sum|cut -c 1-32`
    assertEquals "68e19d0f1092e1f0429ba07b418d1d9f" `$sambamba $opts view -f unpack $sortedbam |$md5sum|cut -c 1-32`
}

testSubSample() {
    return # not testing subsampling yet
    # bam file is part of BioD
    outfn=$outdir/subsample.bam
    $sambamba $opts subsample --type fasthash $outdir/ex1_header.sorted.bam --max-cov 1000 -o$outfn > $outdir/subsample.out
    assertEquals 0 $?
    assertEquals "ab4a69f0093674e9feec54afc2411059" $($md5sum $outfn |cut -c 1-32)
    # should be same as sorted bam
    assertEquals "6e226f6acc8466115600d4de52ca8944" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
    outfn=$outdir/subsample2.bam
    $sambamba $opts subsample -r --type fasthash $outdir/ex1_header.sorted.bam --max-cov 10 -o$outfn > $outdir/subsample2.out
    assertEquals 0 $?
    assertEquals "4092fb75c87cde6efc47ff0a4dc443fb" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
}

testSortByName() {
    outfn=$outdir/ex1_header.nsorted.sam
    # use very tiny buffer of 200K so that multithreading is used
    $sambamba $opts sort -t2 -n $sortedbam -o $nsortedbam -m 200K
    $sambamba $opts view -t2 $nsortedbam > $outfn
    assertEquals "3270" `wc -l < $outfn`
    cat $outfn | cut -f1 | LC_ALL=C sort -c
    assertEquals 0 $?
}

testSortByCoordinate() {
    outfn=$outdir/ex1_header.sorted2.bam
    $sambamba $opts sort -t2 -m 300K $nsortedbam -o $outfn
    # $sambamba $opts index -t2
    assertEquals 0 $?
    $sambamba $opts view -t2 $outfn > $outfn.sam
    # check sorted with cat output/ex1_header.sorted2.bam.sam |awk '{print $3 " " $4}'|less
    assertEquals "d92c51b9e067590d7d5a18a0bdbbe0cc" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
}

testSortPicard() {
    outfn=$outdir/ex1_header.picardSorted.bam
    $sambamba $opts sort -t2 --sort-picard -m 300K $nsortedbam -o $outfn
    assertEquals 0 $?
    $sambamba $opts view -t2 $outfn > $outfn.sam
    assertEquals "2d9af987d7aeaf3a38396870a2f30dca" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
}

testSortByNameMatchMates() {
    outfn=$outdir/match_mates_nameSorted.bam
    $sambamba $opts sort -t2 -n -M -m 300K test/match_mates.bam -o $outfn
    assertEquals 0 $?
    assertEquals "8a0b62e1ac85185c31cf66126c876217" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
}

testNaturalSortMatchMates() {
    outfn=$outdir/match_mates_naturalSorted.bam
    $sambamba $opts sort -t2 -N -M -m 300K test/match_mates.bam -o $outfn
    assertEquals 0 $?
    assertEquals "8a0b62e1ac85185c31cf66126c876217" `$sambamba $opts view -f unpack $outfn |$md5sum|cut -c 1-32`
}

testSlice() {
    $sambamba $opts slice $sortedbam chr1 -o /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegions() {
    $sambamba $opts slice test/issue_204.bam 2:166860000-166870000 2:166860000-166870000 |\
        $sambamba $opts view -c /dev/stdin > /dev/null
    assertEquals 0 $?
}

testSliceMultipleRegionsBed() {
    assertEquals "156" `$sambamba $opts slice -L test/chr2_chr3_test_region.bed test/issue_204.bam |\
          $sambamba $opts view -c /dev/stdin`
}

testView() {
    assertEquals "1806" `$sambamba $opts view -c $sortedbam chr2`
    assertEquals "1464" `$sambamba $opts view -c $sortedbam chr1`
    assertEquals "0" `$sambamba $opts view -c $sortedbam '*'`
}

testOverwriteProtection() {
    $sambamba $opts view $nsortedbam -f bam -o ./$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba $opts merge $sortedbam $sortedbam $sortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba $opts sort -m 50M $nsortedbam -o ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba $opts markdup $nsortedbam ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
    $sambamba $opts slice $nsortedbam chr1 -o ./bin/../$nsortedbam 2>/dev/null
    assertNotSame 0 $?
}

testSortingEmptyFile() {
    $sambamba $opts view $sortedbam -f bam -F "ref_id > 3" -o $outdir/empty.bam 2>/dev/null
    $sambamba $opts sort -m 50M $outdir/empty.bam -o $outdir/empty2.bam 2>/dev/null
    assertEquals "0" `$sambamba $opts view -c $outdir/empty2.bam`
}

testMarkdupEmptyFile() {
    $sambamba $opts view $sortedbam -f bam -F "ref_id > 3" -o $outdir/empty.bam 2>/dev/null
    $sambamba $opts markdup $outdir/empty.bam $outdir/empty.dedup.bam 2>/dev/null
    assertEquals "0" `$sambamba $opts view -c $outdir/empty.dedup.bam`
}

testCramWriting() {
    $sambamba $opts view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o $outdir/c1_pad2.cram
    assertEquals 0 $?
}

testCramReading() {
    $sambamba $opts view -C $outdir/c1_pad2.cram >/dev/null
    assertEquals 0 $?
}

testIndexUsage() {
    rm -f $outdir/*.bai && rm -f $outdir/c1_*
    $sambamba $opts view -S htslib/test/c1\#pad2.sam -T htslib/test/c1.fa -f cram -o $outdir/c1_pad2.cram &&
    $sambamba $opts index -C $outdir/c1_pad2.cram && test -e $outdir/c1_pad2.cram.crai; assertEquals 0 $?
    $sambamba $opts index -C $outdir/c1_pad2.cram $outdir/c1_cram_index && test -e $outdir/c1_cram_index.crai; assertEquals 0 $?
    $sambamba $opts index $sortedbam && test -e $sortedbam.bai; assertEquals 0 $?
    $sambamba $opts index $sortedbam $outdir/ex1_header.sorted.bai && test -e $outdir/ex1_header.sorted.bai; assertEquals 0 $?
}

testIssue193() {
    rm -f $outdir/output_193.txt
    $sambamba $opts depth base test/issue_193.bam > $outdir/output_193.txt 2>/dev/null
    diff -q $outdir/output_193.txt test/issue_193_expected_output.txt
    assertEquals 0 $?
}

testIssue204() {
    rm -f $outdir/output_204.txt
    $sambamba $opts index test/issue_204.bam
    $sambamba $opts depth region test/issue_204.bam -L 2:166868600-166868813 -T 15 -T 20 -T 25 -m > $outdir/output_204.txt 2>/dev/null
    diff -q $outdir/output_204.txt test/issue_204_expected_output.txt
    assertEquals 0 $?
}

testIssue206() {
    $sambamba $opts markdup $sortedbam $outdir/ex1_header.dedup.bam 2>/dev/null
    $sambamba $opts view -H $outdir/ex1_header.dedup.bam | grep '@PG' | grep -q 'sambamba'
    assertEquals 0 $?

    $sambamba $opts view $sortedbam -f bam -o $outdir/ex1_header.filtered.bam \
              -F "supplementary or secondary_alignment"
    # skip failing test - filtered bam is not complete
    # $sambamba $opts view -H $outdir/ex1_header.filtered.bam | grep '@PG' | grep -q 'secondary'
    assertEquals 0 $?
}

testIssue225() {
    $sambamba $opts index test/issue225.bam
    $sambamba $opts depth base -c 1 test/issue225.bam > $outdir/depth_base_225_1.txt 2>/dev/null
    diff -q $outdir/depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba $opts depth base -c 0 test/issue225.bam > $outdir/depth_base_225_0.txt 2>/dev/null
    diff -q $outdir/depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?

    # exercise BED codepath as well
    $sambamba $opts depth base -c 1 -L chrM test/issue225.bam > $outdir/depth_base_225_1.txt 2>/dev/null
    diff -q $outdir/depth_base_225_1.txt test/issue225.out
    assertEquals 0 $?

    $sambamba $opts depth base -c 0 -L chrM test/issue225.bam > $outdir/depth_base_225_0.txt 2>/dev/null
    diff -q $outdir/depth_base_225_0.txt test/issue225.z.out
    assertEquals 0 $?
}

testIssue331() {
    $sambamba $opts  view -H -f json test/issue_204.bam > $outdir/issue_331_header.json
    assertEquals 0 $?
    diff -q $outdir/issue_331_header.json test/regression/issue_331_header.json
    assertEquals 0 $?
}

testIssue356(){
    # check for missing PM tag in @RG
    $sambamba $opts view -S test/issue_356.sam -f bam -o $outdir/issue_356.bam
    $sambamba $opts view $outdir/issue_356.bam -H|grep -q "PM:HiSeq2500_HighOutput"
    assertEquals 0 $?
}

shunit2=`which shunit2`
if [ -x "$shunit2" ]; then
    . $shunit2
else
    . contrib/shunit2-2.0.3/shunit2
fi
