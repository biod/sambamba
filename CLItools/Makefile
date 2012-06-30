FILES=../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d ../sam/serialize.d \
	  ../utils/format.d ../alignmentrange.d ../bamoutput.d ../constants.d \
	  ../bgzfcompress.d ../utils/array.d ../utils/value.d ../region.d \
	  ../jsonserialization.d ../samfile.d ../sam/recordparser.d \
	  ../utils/tagstoragebuilder.d ../utils/msgpack.d

all: sambamba-dmd

../sam/recordparser.d:
	cd ../sam && make recordparser.d

DMDFLAGS=-O -release -inline -g
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g -funroll-all-loops -finline-limit=8192

SAMBAMBA_FILES=sambamba/filter.d sambamba/serializer.d sambamba/sambamba.d

sambamba-dmd: $(FILES) $(SAMBAMBA_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_FILES) -ofbuild/sambamba

sambamba-gdc: $(FILES) $(SAMBAMBA_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_FILES) -o build/sambamba
