FILES=bamfile.d chunkinputstream.d bgzfrange.d \
	  samheader.d reference.d alignment.d \
	  tagstorage.d tagvalue.d utils/switchendianness.d \
	  validation/samheader.d validation/alignment.d utils/algo.d \
	  randomaccessmanager.d virtualoffset.d bai/read.d bai/utils/algo.d \
	  bai/bin.d bai/chunk.d utils/range.d utils/memoize.d sam/serialize.d \
	  utils/format.d alignmentrange.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d

all:
	dmd $(LIBFILES) -oflibbam.so -O -release -inline -shared

debug:
	dmd $(LIBFILES) -oflibbam.so -debug -g -shared

scaffolds:
	dmd utils/switchendianness.d samheader.d reference.d alignment.d tagvalue.d tagstorage.d generate_scaffolds.d -ofgenerate_scaffolds -J.
	./generate_scaffolds

unittests: $(TESTFILES)
	dmd $(TESTFILES) -debug -g -unittest -ofrun_unittests
	./run_unittests

unittests-gdc: $(TESTFILES)
	gdc $(TESTFILES) -O0 -funittest -o run_unittests -lpthread -fdebug
	./run_unittests

test: $(FILES) readbam.d
	dmd $(FILES) readbam.d -ofreadbam -O -release -inline -g

test-gdc: $(FILES) readbam.d
	gdc $(FILES) readbam.d -o readbam -O3 -frelease -fno-bounds-check -fno-assert -lpthread -g -funroll-all-loops -finline-limit=2048

clean:
	rm *.o
