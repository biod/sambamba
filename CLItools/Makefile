FILES=../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d ../sam/serialize.d \
	  ../utils/format.d ../alignmentrange.d ../bamoutput.d ../constants.d \
	  ../bgzfcompress.d ../utils/array.d

DMDFLAGS=-O -release -inline
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g

bam2sam: $(FILES) bam2sam/bam2sam.d
	dmd $(DMDFLAGS) $(FILES) bam2sam/bam2sam.d -ofbuild/bam2sam 

bam2sam-gdc: $(FILES) bam2sam/bam2sam.d
	gdc $(GDCFLAGS) $(FILES) bam2sam/bam2sam.d -o build/bam2sam

bam-fetch: $(FILES) bam-fetch/bam_fetch.d
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) bam-fetch/bam_fetch.d -ofbuild/bam-fetch

bam-fetch-gdc: $(FILES) bam-fetch/bam_fetch.d
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) bam-fetch/bam_fetch.d -o build/bam-fetch
