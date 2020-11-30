# This is a minimalistic make file to build sambamba with ldc2 as per instructions on
# https://github.com/biod/sambamba#compiling-sambamba
#
# Targets (64-bit):
#
#   Linux
#   OSX
#
# Typical usage:
#
#   make LIBRARY_PATH=~/opt/ldc2-$ver-linux-x86_64/lib debug|profile|release|static
#
# Static release with optimization (for releases):
#
#   make LIBRARY_PATH=~/opt/ldc2-$ver-linux-x86_64/lib pgo-static
#
# Debug version
#
#   make LIBRARY_PATH=~/opt/ldc2-$ver-linux-x86_64/lib debug
#

D_COMPILER=ldc2

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SYS = OSX
else
  SYS = LINUX
endif

DFLAGS      = -wi -I. -IBioD -g -J.

DLIBS       = $(LIBRARY_PATH)/libphobos2-ldc.a $(LIBRARY_PATH)/libdruntime-ldc.a
DLIBS_DEBUG = $(LIBRARY_PATH)/libphobos2-ldc-debug.a $(LIBRARY_PATH)/libdruntime-ldc-debug.a
LIBS        = lz4/lib/liblz4.a -L-L$(LIBRARY_PATH) -L-lpthread -L-lm -L-lz
LIBS_STATIC = $(LIBRARY_PATH)/libc.a $(DLIBS) $(LIBRARY_PATH)/liblz4.a -L-lz
SRC         = utils/ldc_version_info_.d utils/lz4.d utils/strip_bcf_header.d $(sort $(wildcard BioD/contrib/undead/*.d BioD/contrib/undead/*/*.d)) utils/version_.d $(sort $(wildcard thirdparty/*.d) $(wildcard BioD/bio/*/*.d BioD/bio/*/*/*.d BioD/bio/*/*/*/*.d BioD/bio/*/*/*/*/*.d) $(wildcard sambamba/*.d sambamba/*/*.d sambamba/*/*/*.d))
OBJ         = $(SRC:.d=.o)
OUT         = bin/sambamba-$(shell cat VERSION)

STATIC_LIB_PATH=-Llz4

.PHONY: all debug release static clean test

all: release

debug:                             DFLAGS += -O0 -d-debug -link-debuglib

profile:                           DFLAGS += -fprofile-instr-generate=profile.raw

coverage:                          DFLAGS += -cov

release static pgo-static:         DFLAGS += -O3 -release -enable-inlining -boundscheck=off -L-lz

static:                            DFLAGS += -static -L-Bstatic

pgo-static:                        DFLAGS += -fprofile-instr-use=profile.data

lz4-static: lz4/lib/liblz4.a

lz4/lib/liblz4.a: lz4/lib/lz4.c lz4/lib/lz4hc.c lz4/lib/lz4frame.c lz4/lib/xxhash.c
	cd lz4/lib && gcc -O3 -c lz4.c lz4hc.c lz4frame.c xxhash.c && $(AR) rcs liblz4.a lz4.o lz4hc.o lz4frame.o xxhash.o

utils/ldc_version_info_.d:
	python3 ./gen_ldc_version_info.py $(shell which ldmd2) > utils/ldc_version_info_.d
	cat utils/ldc_version_info_.d

ldc_version_info: utils/ldc_version_info_.d

build-setup: ldc_version_info lz4-static
	mkdir -p bin/

default debug release static: $(OUT)

coverage: debug

profile: debug
	$(OUT) sort /gnu/data/in_raw.bam -p > /dev/null
	ldc-profdata merge -output=profile.data profile.raw
	rm $(OUT) ./bin/sambamba.o # trigger rebuild

default: all

# ---- Compile step
%.o: %.d
	$(D_COMPILER) $(DFLAGS) -c $< -od=$(dir $@)

singleobj: build-setup
	$(info compile single object...)
	$(D_COMPILER) -singleobj $(DFLAGS) -c -of=$(OUT).o $(SRC)

# ---- Link step
$(OUT): singleobj
	$(info linking...)
	$(D_COMPILER) $(DFLAGS) -of=$(OUT) $(OUT).o $(LINK_OBJ) $(LIBS)

test: $(OUT)
	$(OUT) --version
	./run_tests.sh $(OUT)

check: test

debug-strip:
	objcopy --only-keep-debug bin/sambamba sambamba.debug
	objcopy --strip-debug bin/sambamba
	objcopy --add-gnu-debuglink=sambamba.debug bin/sambamba
	mv sambamba.debug bin/

pgo-static: static debug-strip

install:
	install -m 0755 bin/sambamba $(prefix)/bin

clean: clean-d
	rm lz4/lib/*.[oa]
	rm -f profile.data
	rm -f profile.raw

clean-d:
	rm -rf bin/*
	rm -vf utils/ldc_version_info_.d
	rm -f $(OBJ) $(OUT) trace.{def,log}
