# This is a minimalistic make file to build sambamba with ldc2 as per instructions on
# https://github.com/biod/sambamba#compiling-sambamba

D_COMPILER=ldc2
DFLAGS = -wi -I. -IBioD -IundeaD/src -g

DLIBS  = $(LIBRARY_PATH)/libphobos2-ldc.a $(LIBRARY_PATH)/libdruntime-ldc.a
DLIBS_DEBUG = $(LIBRARY_PATH)/libphobos2-ldc-debug.a $(LIBRARY_PATH)/libdruntime-ldc-debug.a
# RPATH  = -L--rpath=$(dir $(realpath $(LIBRARY_PATH)/libz.so)):$(dir $(realpath $(LIBRARY_PATH)/liblz4.so))
LIBS   = htslib/libhts.a lz4/lib/liblz4.a -L-L$(LIBRARY_PATH) -L-lrt -L-lpthread -L-lm
LIBS_STATIC = $(LIBRARY_PATH)/libc.a $(DLIBS) htslib/libhts.a lz4/lib/liblz4.a
SRC    = $(wildcard main.d utils/*.d thirdparty/*.d cram/*.d) $(wildcard undeaD/src/undead/*.d) $(wildcard BioD/bio/*/*.d BioD/bio/*/*/*.d) $(wildcard sambamba/*.d sambamba/*/*.d sambamba/*/*/*.d)
OBJ    = $(SRC:.d=.o) utils/ldc_version_info_.o
OUT    = build/sambamba

STATIC_LIB_PATH=-Lhtslib -Llz4

.PHONY: all debug release static clean test

debug:       DFLAGS += -O0 -d-debug -link-debuglib

release:     DFLAGS += -O3 -release -enable-inlining -boundscheck=off

all: release

lz4-static: lz4/lib/liblz4.a

lz4/lib/liblz4.a: lz4/lib/lz4.c lz4/lib/lz4hc.c lz4/lib/lz4frame.c lz4/lib/xxhash.c
	cd lz4/lib && $(CC) -O3 -c lz4.c lz4hc.c lz4frame.c xxhash.c && $(AR) rcs liblz4.a lz4.o lz4hc.o lz4frame.o xxhash.o

htslib-static:
	cd htslib && $(MAKE)

ldc-version-info:
	./gen_ldc_version_info.py $(shell which ldmd2) > utils/ldc_version_info_.d

utils/ldc_version_info_.o: ldc-version-info
	$(D_COMPILER) $(DFLAGS) -c utils/ldc_version_info_.d -od=$(dir $@)

build-setup: htslib-static lz4-static ldc-version-info
	mkdir -p build/

default debug release: $(OUT)

default: all

# ---- Compile step
%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

# ---- Link step
$(OUT): build-setup $(OBJ)
	$(D_COMPILER) $(DFLAGS) -of=build/sambamba $(OBJ) $(LIBS)

test:
	./run_tests.sh

check: test

debug-strip:
	objcopy --only-keep-debug build/sambamba sambamba.debug
	objcopy --strip-debug build/sambamba
	objcopy --add-gnu-debuglink=sambamba.debug build/sambamba
	mv sambamba.debug build/

install:
	install -m 0755 build/sambamba $(prefix)/bin

clean: clean-d
	cd htslib ; make clean

clean-d:
	rm -rf build/*
	rm -f $(OBJ) $(OUT) trace.{def,log}
