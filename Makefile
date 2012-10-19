FILES=constants.d \
	utils/array.d \
	utils/range.d \
	utils/stream.d \
	utils/memoize.d \
	utils/format.d \
	utils/algo.d \
	utils/msgpack.d \
	utils/switchendianness.d \
	tagvalue.d \
	utils/tagstoragebuilder.d \
	utils/value.d \
	virtualoffset.d \
	bai/chunk.d \
	bai/bin.d \
	bai/read.d \
	bai/utils/algo.d \
	BioD/TinyMap.d \
	BioD/Base.d \
	bgzfblock.d \
	chunkinputstream.d \
	alignment.d \
	alignmentrange.d \
	bgzfrange.d \
	randomaccessmanager.d \
	reference.d \
	samheader.d \
	bamfile.d \
	validation/samheader.d \
	validation/alignment.d \
	utils/graph.d \
	utils/samheadermerger.d \
	sam/serialize.d \
	sam/recordparser.d \
	samfile.d \
	md/operation.d \
	md/parse.d \
	md/core.d \
	reconstruct.d \
	splitter.d \
	pileuprange.d \
	bgzfcompress.d \
	bamoutput.d

LIBFILES = $(FILES) bindings.d
TESTFILES = $(FILES) unittests.d utils/tmpfile.d

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
