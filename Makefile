FILES=bamfile.d rangetransformer.d chunkinputstream.d bgzfrange.d \
	  utils/inputrangechunks.d samheader.d reference.d alignment.d \
	  tagstorage.d tagvalue.d utils/switchendianness.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d

all: scaffolds
	dmd $(LIBFILES) -oflibbam.so -O -release -inline -shared

debug:
	dmd $(LIBFILES) -oflibbam.so -debug -g -shared

scaffolds:
	dmd samheader.d reference.d generate_scaffolds.d -ofgenerate_scaffolds -J.
	./generate_scaffolds

unittests: $(TESTFILES)
	dmd $(TESTFILES) -debug -g -unittest -ofrun_unittests
	./run_unittests

unittests-gdc: $(TESTFILES)
	/opt/gdc/bin/gdc $(TESTFILES) -funittest -o run_unittests -lpthread
	./run_unittests

test: $(FILES) readbam.d
	dmd $(FILES) readbam.d -ofreadbam -O -release -inline

test-gdc: $(FILES) readbam.d
	/opt/gdc/bin/gdc $(FILES) readbam.d -o readbam -O3 -lpthread -frelease -fno-bounds-check -fno-assert -fversion=serial

clean:
	rm *.o
