# This is a minimalistic make file to build sambamba with ldc2 as per instructions on
# https://github.com/biod/sambamba#compiling-sambamba
#
# Typical usage:
#
#    make LIBRARY_PATH=~/opt/ldc2-1.7.0-linux-x86_64/lib debug|profile|release|static
#
# Static release with optimization (for releases):
#
#   make LIBRARY_PATH=~/opt/ldc2-1.7.0-linux-x86_64/lib pgo-static
#
# Debug version
#
#   make LIBRARY_PATH=~/opt/ldc2-1.7.0-linux-x86_64/lib debug

D_COMPILER=ldc2
DFLAGS      = -wi -I. -IBioD -IundeaD/src -g

DLIBS       = $(LIBRARY_PATH)/libphobos2-ldc.a $(LIBRARY_PATH)/libdruntime-ldc.a
DLIBS_DEBUG = $(LIBRARY_PATH)/libphobos2-ldc-debug.a $(LIBRARY_PATH)/libdruntime-ldc-debug.a
LIBS        = htslib/libhts.a lz4/lib/liblz4.a -L-L$(LIBRARY_PATH) -L-lrt -L-lpthread -L-lm
LIBS_STATIC = $(LIBRARY_PATH)/libc.a $(DLIBS) htslib/libhts.a lz4/lib/liblz4.a
SRC         = $(wildcard main.d utils/*.d thirdparty/*.d cram/*.d) $(wildcard undeaD/src/undead/*.d) $(wildcard BioD/bio/*/*.d BioD/bio/*/*/*.d BioD/bio2/*.d BioD/bio2/*/*.d) $(wildcard sambamba/*.d sambamba/*/*.d sambamba/*/*/*.d)
OBJ         = $(SRC:.d=.o) utils/ldc_version_info_.o
OUT         = bin/sambamba

STATIC_LIB_PATH=-Lhtslib -Llz4

.PHONY: all debug release static clean test

debug:                             DFLAGS += -O0 -d-debug -link-debuglib

profile:                           DFLAGS += -fprofile-instr-generate=profile.raw

release static profile pgo-static: DFLAGS += -O3 -release -enable-inlining -boundscheck=off

static:                            DFLAGS += -static -L-Bstatic

pgo-static:                        DFLAGS += -fprofile-instr-use=profile.data

all: release

lz4-static: lz4/lib/liblz4.a

lz4/lib/liblz4.a: lz4/lib/lz4.c lz4/lib/lz4hc.c lz4/lib/lz4frame.c lz4/lib/xxhash.c
	cd lz4/lib && $(CC) -O3 -c lz4.c lz4hc.c lz4frame.c xxhash.c && $(AR) rcs liblz4.a lz4.o lz4hc.o lz4frame.o xxhash.o

htslib-static:
	cd htslib && $(MAKE)

ldc-version-info:
	./gen_ldc_version_info.py $(shell which ldmd2) > utils/ldc_version_info_.d
	cat utils/ldc_version_info_.d

utils/ldc_version_info_.o: ldc-version-info
	$(D_COMPILER) $(DFLAGS) -c utils/ldc_version_info_.d -od=$(dir $@)

build-setup: htslib-static lz4-static ldc-version-info
	mkdir -p bin/

default debug release static: $(OUT)

profile: release
	./bin/sambamba sort /gnu/data/in_raw.bam -p > /dev/null
	ldc-profdata merge -output=profile.data profile.raw
	rm ./bin/sambamba ./bin/sambamba.o # trigger rebuild

default: all

# ---- Compile step
%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

singleobj:
	$(info compile single object...)
	$(D_COMPILER) -singleobj $(DFLAGS) -c -of=bin/sambamba.o $(SRC)

# ---- Link step
$(OUT): build-setup singleobj utils/ldc_version_info_.o
	$(info linking...)
	$(D_COMPILER) $(DFLAGS) -of=bin/sambamba bin/sambamba.o utils/ldc_version_info_.o $(LIBS)

test:
	./run_tests.sh

check: test

debug-strip:
	objcopy --only-keep-debug bin/sambamba sambamba.debug
	objcopy --strip-debug bin/sambamba
	objcopy --add-gnu-debuglink=sambamba.debug bin/sambamba
	mv sambamba.debug bin/

pgo-static: profile static debug-strip

install:
	install -m 0755 bin/sambamba $(prefix)/bin

clean: clean-d
	cd htslib ; make clean
	rm -f profile.data
	rm -f profile.raw

clean-d:
	rm -rf bin/*
	rm -f $(OBJ) $(OUT) trace.{def,log}
