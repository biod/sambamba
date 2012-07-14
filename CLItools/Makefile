FILES=../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d ../sam/serialize.d \
	  ../utils/format.d ../alignmentrange.d ../bamoutput.d ../constants.d \
	  ../bgzfcompress.d ../utils/array.d ../utils/value.d ../region.d \
	  ../jsonserialization.d ../bai/indexing.d

SAMFILES=../samfile.d ../sam/recordparser.d ../utils/tagstoragebuilder.d

all: sambamba-dmd

../sam/recordparser.d:
	cd ../sam && make recordparser.d

DMDFLAGS=-O -release -inline -g
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g -funroll-all-loops -finline-limit=8192

SAMBAMBA_FILES=sambamba/filter.d sambamba/serializer.d sambamba/sambamba.d
SAMBAMBA_INDEX_FILES=sambamba-index/index.d
SAMBAMBA_SORT_MERGE_COMMON_FILES=common/comparators.d common/nwayunion.d
SAMBAMBA_SORT_FILES=sambamba-sort/sort.d sambamba-sort/thirdparty/mergesort.d $(SAMBAMBA_SORT_MERGE_COMMON_FILES)
SAMBAMBA_MERGE_FILES=sambamba-merge/merge.d ../utils/samheadermerger.d ../utils/graph.d $(SAMBAMBA_SORT_MERGE_COMMON_FILES)
SAMBAMBA_FLAGSTAT_FILES=sambamba-flagstat/flagstat.d

sambamba-dmd: $(FILES) $(SAMFILES) $(SAMBAMBA_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMFILES) $(SAMBAMBA_FILES) -ofbuild/sambamba

sambamba-gdc: $(FILES) $(SAMFILES) $(SAMBAMBA_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMFILES) $(SAMBAMBA_FILES) -o build/sambamba

sambamba-index-dmd: $(FILES) $(SAMBAMBA_INDEX_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_INDEX_FILES) -ofbuild/sambamba-index

sambamba-index-gdc: $(FILES) $(SAMBAMBA_INDEX_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_INDEX_FILES) -o build/sambamba-index

sambamba-sort-dmd: $(FILES) $(SAMBAMBA_SORT_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_SORT_FILES) -ofbuild/sambamba-sort

sambamba-sort-gdc: $(FILES) $(SAMBAMBA_SORT_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_SORT_FILES) -o build/sambamba-sort

sambamba-merge-dmd: $(FILES) $(SAMBAMBA_MERGE_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_MERGE_FILES) -ofbuild/sambamba-merge

sambamba-merge-gdc: $(FILES) $(SAMBAMBA_MERGE_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_MERGE_FILES) -o build/sambamba-merge

sambamba-flagstat-dmd: $(FILES) $(SAMBAMBA_FLAGSTAT_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_FLAGSTAT_FILES) -ofbuild/sambamba-flagstat

sambamba-flagstat-gdc: $(FILES) $(SAMBAMBA_FLAGSTAT_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_FLAGSTAT_FILES) -o build/sambamba-flagstat

sambamba-dmd-dev: $(FILES) $(SAMFILES) $(SAMBAMBA_FILES)
	mkdir -p build
	dmd -g $(FILES) $(SAMFILES) $(SAMBAMBA_FILES) -version=development -ofbuild/sambamba

sambamba-gdc-dev: $(FILES) $(SAMFILES) $(SAMBAMBA_FILES)
	mkdir -p build
	gdc -g -O0 -lpthread $(FILES) $(SAMFILES) $(SAMBAMBA_FILES) -fversion=development -o build/sambamba
