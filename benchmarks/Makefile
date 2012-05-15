FILES=../bamfile.d ../rangetransformer.d ../chunkinputstream.d ../bgzfrange.d \
	  ../utils/inputrangechunks.d ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d

serial: $(FILES) benchmarks.d
	dmd $(FILES) benchmarks.d -ofrun_benchmarks -O -release -inline -version=serial -version=measure

parallel: $(FILES) benchmarks.d
	dmd $(FILES) benchmarks.d -ofrun_benchmarks_mt -O -release -inline -version=parallel -version=measure

clean:
	rm *.o
