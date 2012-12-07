#!/usr/bin/env bash

rm -rf build/
mkdir -p build/
cd build/

FILES="../main.d ../BioD/bio/sam/utils/fastrecordparser.d `rdmd --makedepend -I../BioD/ ../main.d | cut -d : -f 2 | awk -v RS=' ' -v ORS=' ' '!/dmd/ {print}'`"

GDC_RELEASE_FLAGS='-g -O3 -finline-limit=8192 -funroll-all-loops -frelease -finline -fno-bounds-check -fno-assert'
GDC_DEBUG_FLAGS='-g'

gdc -I.. -I../BioD -m$1 $GDC_RELEASE_FLAGS -lpthread $FILES -o sambamba
