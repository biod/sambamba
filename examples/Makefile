FILES=../bamfile.d ../chunkinputstream.d ../bgzfrange.d \
	  ../samheader.d ../reference.d ../alignment.d \
	  ../tagstorage.d ../tagvalue.d ../utils/switchendianness.d \
	  ../validation/alignment.d ../utils/algo.d ../randomaccessmanager.d \
	  ../bai/read.d ../bai/chunk.d ../bai/bin.d ../bai/utils/algo.d \
	  ../virtualoffset.d ../utils/range.d ../utils/memoize.d ../sam/serialize.d \
	  ../utils/format.d ../alignmentrange.d ../constants.d ../bamoutput.d \
	  ../bgzfcompress.d

DMDFLAGS=-O -release -inline
GDCFLAGS=-O3 -frelease -finline -fno-assert -fno-bounds-check -lpthread -g

validator: $(FILES) validator.d
	dmd $(FILES) validator.d -ofvalidator $(DMDFLAGS) -version=serial

validator-gdc: $(FILES) validator.d
	gdc $(FILES) validator.d -o validator $(GDCFLAGS) -fversion=serial

randomaccess: $(FILES) randomaccess.d
	dmd $(FILES) randomaccess.d -ofrandomaccess $(DMDFLAGS) -version=serial

randomaccess-gdc: $(FILES) randomaccess.d
	gdc $(FILES) randomaccess.d -o randomaccess $(GDCFLAGS) -fversion=serial

clean:
	rm *.o
