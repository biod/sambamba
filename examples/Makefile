FILES=../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d

DMDFLAGS=-O -release -inline -g
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g

bam2sam: $(FILES) bam2sam.d
	dmd $(FILES) bam2sam.d -ofbam2sam $(DMDFLAGS) -version=serial

bam2sam-gdc: $(FILES) bam2sam.d
	gdc $(FILES) bam2sam.d -o bam2sam $(GDCFLAGS) -fversion=serial

validator: $(FILES) validator.d
	dmd $(FILES) validator.d -ofvalidator $(DMDFLAGS) -version=serial

validator-gdc: $(FILES) validator.d
	gdc $(FILES) validator.d -o validator $(GDCFLAGS) -fversion=serial

randomaccess: $(FILES) randomaccess.d
	dmd $(FILES) randomaccess.d -ofrandomaccess $(DMDFLAGS) -version=serial

randomaccess-gdc: $(FILES) randomaccess.d
	gdc $(FILES) randomaccess.d -o randomaccess $(GDCFLAGS) -fversion=serial

clean:
	rm *.o
