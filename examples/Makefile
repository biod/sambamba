FILES=../bamfile.d ../rangetransformer.d ../chunkinputstream.d ../bgzfrange.d \
	  ../utils/inputrangechunks.d ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d

bam2sam: $(FILES) bam2sam.d
	dmd $(FILES) bam2sam.d -ofbam2sam -O -release -inline -version=serial -g

bam2sam-gdc: $(FILES) bam2sam.d
	gdc $(FILES) bam2sam.d -o bam2sam -O3 -frelease -finline -fno-bounds-check -fno-assert -fversion=parallel -g -lpthread

validator: $(FILES) validator.d
	dmd $(FILES) validator.d -ofvalidator -O -inline -version=serial -g

validator-gdc: $(FILES) validator.d
	gdc $(FILES) validator.d -ofvalidator -O3 -finline -frelease -fno-bounds-check -fno-assert -g -lpthread

clean:
	rm *.o
