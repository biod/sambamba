FILES=bamfile.d chunkinputstream.d bgzfrange.d \
	  samheader.d reference.d alignment.d \
	  tagvalue.d utils/switchendianness.d \
	  validation/samheader.d validation/alignment.d utils/algo.d \
	  randomaccessmanager.d virtualoffset.d bai/read.d bai/utils/algo.d \
	  bai/bin.d bai/chunk.d utils/range.d utils/memoize.d sam/serialize.d \
	  utils/format.d alignmentrange.d bamoutput.d constants.d bgzfcompress.d \
	  utils/array.d utils/value.d samfile.d sam/recordparser.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d

FILESTODOCUMENT = bamfile.d alignment.d reference.d tagvalue.d \
				  samheader.d validation/samheader.d validation/samheader.d

all:
	dmd $(LIBFILES) -oflibbam.so -O -release -inline -shared

debug:
	dmd $(LIBFILES) -oflibbam.so -debug -g -shared

scaffolds:
	dmd utils/switchendianness.d samheader.d reference.d alignment.d tagvalue.d tagstorage.d generate_scaffolds.d -ofgenerate_scaffolds -J.
	./generate_scaffolds

region-parser: region.rl
	ragel region.rl -D -G2

unittests: $(TESTFILES)
	dmd $(TESTFILES) -debug -g -unittest -ofrun_unittests -version=serial
	./run_unittests

unittests-gdc: $(TESTFILES)
	gdc $(TESTFILES) -O0 -funittest -o run_unittests -lpthread -fdebug -fversion=serial
	./run_unittests

test: $(FILES) readbam.d
	dmd $(FILES) readbam.d -ofreadbam -g -debug

test-gdc: $(FILES) readbam.d
	gdc $(FILES) readbam.d -o readbam -g -fdebug -lpthread

testsam: $(FILES) readsam.d samfile.d sam/sam_alignment.d
	rdmd --compiler=gdmd -O -release -inline --build-only -L-lpthread -g readsam.d

clean:
	rm *.o
