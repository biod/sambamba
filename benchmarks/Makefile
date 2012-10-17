FILES=../bamfile.d ../utils/range.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagvalue.d ../utils/switchendianness.d \
	  ../virtualoffset.d ../randomaccessmanager.d ../bai/read.d \
	  ../bai/bin.d ../bai/chunk.d ../bai/utils/algo.d ../utils/algo.d \
	  ../alignmentrange.d ../utils/memoize.d ../constants.d ../utils/array.d \
	  ../utils/value.d ../bgzfcompress.d ../utils/format.d ../BioD/Base.d \
	  ../BioD/TinyMap.d ../utils/msgpack.d ../utils/stream.d


DMDFLAGS=-O -release -inline
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread

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
