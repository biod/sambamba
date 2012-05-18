FILES=../bamfile.d ../rangetransformer.d ../chunkinputstream.d ../bgzfrange.d \
	  ../utils/inputrangechunks.d ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d

DMDFLAGS=-O -release -inline
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check

all: serial parallel serial-gdc parallel-gdc

serial: $(FILES) benchmarks.d
	dmd $(FILES) benchmarks.d -ofrun_benchmarks $(DMDFLAGS) -version=serial

parallel: $(FILES) benchmarks.d
	dmd $(FILES) benchmarks.d -ofrun_benchmarks_mt $(DMDFLAGS) -version=parallel

serial-gdc: $(FILES) benchmarks.d
	gdc $(FILES) benchmarks.d -o run_benchmarks_gdc $(GDCFLAGS) -fversion=serial 

parallel-gdc: $(FILES) benchmarks.d
	gdc $(FILES) benchmarks.d -o run_benchmarks_mt_gdc $(GDCFLAGS) -fversion=parallel

clean:
	rm *.o
