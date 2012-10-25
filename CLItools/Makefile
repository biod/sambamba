FILES=$(addprefix ../,\
	constants.d \
	utils/array.d \
	utils/range.d \
	utils/stream.d \
	utils/memoize.d \
	utils/format.d \
	utils/algo.d \
	utils/msgpack.d \
	utils/switchendianness.d \
	tagvalue.d \
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
	sam/serialize.d \
	md/operation.d \
	md/parse.d \
	md/core.d \
	reconstruct.d \
	bgzfcompress.d \
	bamoutput.d \
	bai/indexing.d \
	jsonserialization.d \
	region.d \
	CLItools/common/progressbar.d \
	splitter.d \
	validation/alignment.d)

SAMFILES=../samfile.d ../sam/recordparser.d ../utils/tagstoragebuilder.d

all: sambamba-dmd

../sam/recordparser.d:
	cd ../sam && make recordparser.d

DMDFLAGSFULL=-g -O -release -inline#-O -release -inline -g
GDCFLAGSFULL=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g -funroll-all-loops -finline-limit=8192

DMDFLAGS=$(DMDFLAGSFULL) -version=standalone
GDCFLAGS=$(GDCFLAGSFULL) -fversion=standalone

SAMBAMBA_VIEW_FILES=sambamba-view/filtering.d sambamba-view/headerserializer.d sambamba-view/alignmentrangeprocessor.d sambamba-view/view.d sambamba-view/queryparser.d ../utils/pratt_parser.d
SAMBAMBA_INDEX_FILES=sambamba-index/index.d
SAMBAMBA_SORT_MERGE_COMMON_FILES=common/comparators.d common/nwayunion.d
SAMBAMBA_SORT_FILES=sambamba-sort/sort.d sambamba-sort/thirdparty/mergesort.d ../utils/tmpfile.d
SAMBAMBA_MERGE_FILES=sambamba-merge/merge.d ../utils/samheadermerger.d ../utils/graph.d ../validation/samheader.d
SAMBAMBA_FLAGSTAT_FILES=sambamba-flagstat/flagstat.d

SAMBAMBA_FILES=$(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES) $(SAMBAMBA_INDEX_FILES) $(SAMBAMBA_FLAGSTAT_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) $(SAMBAMBA_SORT_FILES) $(SAMBAMBA_MERGE_FILES) sambamba/main.d

sambamba-dmd: $(SAMBAMBA_FILES)
	mkdir -p build
	dmd $(DMDFLAGSFULL) $(SAMBAMBA_FILES) -ofbuild/sambamba

sambamba-gdc: $(SAMBAMBA_FILES) fastrecordparser
	mkdir -p build
	gdc $(GDCFLAGSFULL) $(SAMBAMBA_FILES) -o build/sambamba

sambamba-gdc-32: $(SAMBAMBA_FILES) fastrecordparser
	mkdir -p build
	gdc -m32 $(GDCFLAGSFULL) $(SAMBAMBA_FILES) -o build/sambamba

sambamba-dmd-serial: $(SAMBAMBA_FILES)
	mkdir -p build
	dmd $(DMDFLAGSFULL) $(SAMBAMBA_FILES) -version=serial -ofbuild/sambamba

sambamba-view-dmd: $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES) -ofbuild/sambamba-view

sambamba-view-gdc: $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES) -o build/sambamba-view

sambamba-index-dmd: $(FILES) $(SAMBAMBA_INDEX_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_INDEX_FILES) -ofbuild/sambamba-index

sambamba-index-gdc: $(FILES) $(SAMBAMBA_INDEX_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_INDEX_FILES) -o build/sambamba-index

sambamba-sort-dmd: $(FILES) $(SAMBAMBA_SORT_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_SORT_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) -ofbuild/sambamba-sort

sambamba-sort-gdc: $(FILES) $(SAMBAMBA_SORT_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_SORT_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) -o build/sambamba-sort

sambamba-merge-dmd: $(FILES) $(SAMBAMBA_MERGE_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) 
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_MERGE_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) -ofbuild/sambamba-merge

sambamba-merge-gdc: $(FILES) $(SAMBAMBA_MERGE_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES)  
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_MERGE_FILES) $(SAMBAMBA_SORT_MERGE_COMMON_FILES) -o build/sambamba-merge

sambamba-flagstat-dmd: $(FILES) $(SAMBAMBA_FLAGSTAT_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMBAMBA_FLAGSTAT_FILES) -ofbuild/sambamba-flagstat

sambamba-flagstat-gdc: $(FILES) $(SAMBAMBA_FLAGSTAT_FILES)
	mkdir -p build
	gdc $(GDCFLAGS) $(FILES) $(SAMBAMBA_FLAGSTAT_FILES) -o build/sambamba-flagstat

sambamba-dmd-dev: $(SAMBAMBA_FILES)
	mkdir -p build
	dmd -g $(SAMBAMBA_FILES) -version=development -ofbuild/sambamba

sambamba-gdc-dev: $(SAMBAMBA_FILES)
	mkdir -p build
	gdc -g -O0 -lpthread $(SAMBAMBA_FILES) -fversion=development -o build/sambamba

fastrecordparser: ../sam/sam_alignment.rl
	cd ../sam && make fastrecordparser && cd ..

sambamba-view-dmd-serial: $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES)
	mkdir -p build
	dmd $(DMDFLAGS) $(FILES) $(SAMFILES) $(SAMBAMBA_VIEW_FILES) -version=serial -ofbuild/sambamba-view

sambamba-fetchregion-dmd: $(FILES) sambamba-fetchregion/fetchregion.d
	rdmd --force --build-only -ofbuild/sambamba-fetchregion -O -release -inline -I.. sambamba-fetchregion/fetchregion.d

clean:
	rm -rf build/
