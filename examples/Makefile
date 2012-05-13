FILES=../bamfile.d ../rangetransformer.d ../chunkinputstream.d ../bgzfrange.d \
	  ../utils/inputrangechunks.d ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d

bam2sam: $(FILES) bam2sam.d
	dmd $(FILES) bam2sam.d -ofbam2sam -O -release -inline -version=serial

clean:
	rm *.o
