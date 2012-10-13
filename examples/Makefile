FILES=../utils/stream.d ../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../BioD/Base.d ../alignment.d \
	  ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d ../sam/serialize.d \
	  ../utils/format.d ../alignmentrange.d ../constants.d ../bamoutput.d \
	  ../bgzfcompress.d ../utils/array.d ../utils/value.d ../utils/msgpack.d \
	  ../reconstruct.d ../md/core.d ../md/operation.d ../md/parse.d

DMDFLAGS=-g -O -release -inline#-g -debug
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g -funroll-all-loops -finline-limit=8192 -I..

DCOLLECTIONS_FILES=../dcollections/LinkList.d ../dcollections/model/List.d \
				   ../dcollections/model/Addable.d ../dcollections/model/Iterator.d \
				   ../dcollections/Link.d ../dcollections/DefaultFunctions.d \
				   ../dcollections/DefaultAllocator.d

coverage: $(FILES) coverage.d ../pileuprange.d
	dmd $(FILES) coverage.d ../pileuprange.d $(DCOLLECTIONS_FILES) $(DMDFLAGS) -ofcoverage

# segfaults :-((( perhaps, it's possible to tweak the file order somehow...
coverage-gdc: $(FILES) coverage.d ../pileuprange.d
	gdc $(FILES) $(DCOLLECTIONS_FILES) ../pileuprange.d coverage.d $(GDCFLAGS) -o coverage

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
