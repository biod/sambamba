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
	bgzfcompress.d \
	bamoutput.d \
	splitter.d)

DMDFLAGS=-g -O -release -inline#-g -debug
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g -funroll-all-loops -finline-limit=8192 -I..

coverage: $(FILES) coverage.d ../pileuprange.d
	dmd $(FILES) coverage.d ../pileuprange.d $(DMDFLAGS) -ofcoverage

coverage-gdc: $(FILES) coverage.d ../pileuprange.d
	gdc $(FILES) ../pileuprange.d coverage.d $(GDCFLAGS) -o coverage

snpcaller: $(FILES) ../pileuprange.d ../BioD/Call.d ../BioD/TinyMap.d ../BioD/Genotype.d ../snpcallers/maq.d snpcaller.d
	rdmd --build-only  -I.. $(DMDFLAGS) snpcaller.d

validator: $(FILES) validator.d
	dmd $(FILES) validator.d -ofvalidator $(DMDFLAGS) -version=serial

validator-gdc: $(FILES) validator.d
	gdc $(FILES) validator.d -o validator $(GDCFLAGS)

randomaccess: $(FILES) randomaccess.d
	dmd $(FILES) randomaccess.d -ofrandomaccess $(DMDFLAGS) -version=serial

randomaccess-gdc: $(FILES) randomaccess.d
	gdc $(FILES) randomaccess.d -o randomaccess $(GDCFLAGS) -fversion=serial

rename: $(FILES) rename.d
	dmd $(FILES) rename.d -ofrenamereads $(DMDFLAGS) -version=parallel

rename-gdc: $(FILES) rename.d
	gdc $(FILES) rename.d -o renamereads $(GDCFLAGS) -fversion=serial

printindex: $(FILES) printindex.d
	dmd $(FILES) printindex.d -ofprintindex $(DMDFLAGS) -version=serial

printindex-gdc: $(FILES) printindex.d
	gdc $(FILES) printindex.d -o printindex $(GDCFLAGS) -fversion=serial

msgpackoutput: $(FILES) msgpackoutput.d 
	dmd $(FILES) msgpackoutput.d -ofmsgpackoutput $(DMDFLAGS)

msgpackoutput-gdc: $(FILES) msgpackoutput.d
	gdc $(FILES) msgpackoutput.d -o msgpackoutput $(GDCFLAGS)

clean:
	rm *.o
