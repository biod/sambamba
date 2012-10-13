FILES=utils/stream.d bamfile.d chunkinputstream.d bgzfrange.d \
	  samheader.d reference.d BioD/Base.d alignment.d \
	  tagvalue.d utils/switchendianness.d \
	  validation/samheader.d validation/alignment.d utils/algo.d \
	  randomaccessmanager.d virtualoffset.d bai/read.d bai/utils/algo.d \
	  bai/bin.d bai/chunk.d utils/range.d utils/memoize.d sam/serialize.d \
	  utils/format.d alignmentrange.d bamoutput.d constants.d bgzfcompress.d \
	  utils/array.d utils/value.d utils/tagstoragebuilder.d samfile.d \
	  sam/recordparser.d utils/samheadermerger.d utils/graph.d utils/msgpack.d \
	  reconstruct.d md/core.d md/operation.d md/parse.d splitter.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d utils/tmpfile.d

FILESTODOCUMENT = bamfile.d alignment.d reference.d tagvalue.d \
				  samheader.d validation/samheader.d validation/samheader.d

all: unittests

sam/recordparser.d : sam/sam_alignment.rl
	cd sam && make recordparser.d

debug:
	dmd $(LIBFILES) -oflibbam.so -debug -g -shared

region-parser: region.rl
	ragel region.rl -D -G2

unittests: $(TESTFILES)
	dmd $(TESTFILES) -debug -g -unittest -ofrun_unittests -version=serial
	./run_unittests

unittests-gdc: $(TESTFILES)
	gdc $(TESTFILES) -O0 -g -funittest -o run_unittests -lpthread -fdebug -fversion=serial
	./run_unittests

clean:
	rm -f *.o 
	cd sam && make clean
