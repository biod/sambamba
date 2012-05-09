FILES=bamfile.d rangetransformer.d chunkinputstream.d bgzfrange.d \
	  utils/inputrangechunks.d samheader.d reference.d alignment.d

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

test: $(FILES) test_tags.d
	dmd $(FILES) test_tags.d -ofreadbam -O -release -inline

clean:
	rm *.o
