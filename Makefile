FILES=bamfile.d rangetransformer.d chunkinputstream.d bgzfrange.d \
	  utils/inputrangechunks.d samheader.d reference.d alignment.d \
	  tagstorage.d tagvalue.d utils/switchendianness.d \
	  validation/samheader.d utils/algo.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d

all: scaffolds
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
	gdc $(TESTFILES) -funittest -o run_unittests -lpthread
	./run_unittests

test: $(FILES) readbam.d
	dmd $(FILES) readbam.d -ofreadbam -O -release -inline -version=serial -g

test-gdc: $(FILES) readbam.d
	gdc $(FILES) readbam.d -o readbam -O3 -frelease -fno-bounds-check -fno-assert -lpthread -fversion=serial -g

clean:
	rm *.o
