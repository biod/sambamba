D_COMPILER=dmd
D_FLAGS=--compiler=dmd -IBioD -g -d#-O -release -inline # -version=serial
LDMD=ldmd2

STATIC_LIB_PATH=-Lhtslib -Llz4/lib
STATIC_LIB_SUBCMD=$(STATIC_LIB_PATH) -Wl,-Bstatic -lhts -llz4 -Wl,-Bdynamic
RDMD_FLAGS=--force --build-only --compiler=$(D_COMPILER) $(D_FLAGS)

PLATFORM := $(shell uname -s)

ifeq "$(PLATFORM)" "Darwin"

LINK_CMD=gcc -dead_strip -lphobos2-ldc -ldruntime-ldc -lm -lpthread htslib/libhts.a lz4/lib/liblz4.a build/sambamba.o -o build/sambamba
DMD_STATIC_LIBS=htslib/libhts.a lz4/lib/liblz4.a

define split-debug
dsymutil build/sambamba -o build/sambamba.dSYM
strip -S build/sambamba
endef

else

LINK_CMD=gcc -Wl,--gc-sections -o build/sambamba build/sambamba.o $(STATIC_LIB_SUBCMD) -l:libphobos2-ldc.a -l:libdruntime-ldc.a  -lrt -lpthread -lm
DMD_STATIC_LIBS=-L-Lhtslib -L-l:libhts.a -L-l:libphobos2.a -L-Llz4/lib -L-l:liblz4.a

define split-debug
objcopy --only-keep-debug build/sambamba sambamba.debug
objcopy --strip-debug build/sambamba
objcopy --add-gnu-debuglink=sambamba.debug build/sambamba
mv sambamba.debug build/
endef

endif

PREREQS := ldc-version-info htslib-static lz4-static

# DMD only - this goal is used because of fast compilation speed, during development
all: $(PREREQS)
	mkdir -p build/
	rdmd --force --build-only $(D_FLAGS) $(DMD_STATIC_LIBS) -ofbuild/sambamba main.d

# This is the main Makefile goal, used for building releases (best performance)
sambamba-ldmd2-64: $(PREREQS)
	mkdir -p build/
	$(LDMD) @sambamba-ldmd-release.rsp
	$(LINK_CMD)
	$(split-debug)

# For debugging; GDB & Valgrind are more friendly to executables created using LDC/GDC than DMD
sambamba-ldmd2-debug: $(PREREQS)
	mkdir -p build/
	$(LDMD) @sambamba-ldmd-debug.rsp
	$(LINK_CMD)

ldc-version-info:
	./gen_ldc_version_info.py $(shell which $(LDMD)) > utils/ldc_version_info_.d

htslib-static:
	cd htslib && $(MAKE)

lz4-static: lz4/lib/liblz4.a

lz4/lib/liblz4.a: lz4/lib/lz4.c lz4/lib/lz4hc.c lz4/lib/lz4frame.c lz4/lib/xxhash.c
	cd lz4/lib && $(CC) -O3 -c lz4.c lz4hc.c lz4frame.c xxhash.c && $(AR) rcs liblz4.a lz4.o lz4hc.o lz4frame.o xxhash.o

# all below link to libhts dynamically for simplicity

sambamba-flagstat:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-flagstat sambamba/flagstat.d

sambamba-merge:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-merge sambamba/merge.d

sambamba-index:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-index sambamba/index.d

sambamba-sort:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-sort sambamba/sort.d

sambamba-view:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-view sambamba/view.d

sambamba-slice:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-slice sambamba/slice.d

sambamba-markdup:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-markdup sambamba/markdup.d

sambamba-depth:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-depth sambamba/depth.d

sambamba-pileup:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -L-lhts -version=standalone -ofbuild/sambamba-pileup sambamba/pileup.d

.PHONY: clean ldc-version-info

clean:
	rm -rf build/ ; $(MAKE) -C htslib clean ; $(MAKE) -C lz4 clean
