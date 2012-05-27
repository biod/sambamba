FILES=../bamfile.d ../rangetransformer.d ../chunkinputstream.d ../bgzfrange.d \
	  ../utils/inputrangechunks.d ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d

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

voffset: $(FILES) voffset.d
	dmd $(FILES) voffset.d -ofvoffset -version=serial -g

clean:
	rm *.o
