D_COMPILER=ldmd2
D_FLAGS=-IBioD -O2 -release #-inline -g# -version=serial

RDMD_FLAGS=--force --compiler=$(D_COMPILER) --build-only $(D_FLAGS)

all:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -ofbuild/sambamba main.d

release:
	mkdir -p build/
	rdmd --build-only --force -O -release -inline -IBioD/ -ofbuild/sambamba main.d

sambamba-ldmd2-64:
	mkdir -p build/
	rdmd --build-only --force --compiler=ldmd2 -O2 -m64 -noboundscheck -release -inline -IBioD/ -ofbuild/sambamba main.d

sambamba-ldmd2-32:
	mkdir -p build/
	rdmd --build-only --force --compiler=ldmd2 -O2 -m32 -release -inline -IBioD/ -ofbuild/sambamba main.d

sambamba-flagstat:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-flagstat sambamba/flagstat.d

sambamba-merge:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-merge sambamba/merge.d

sambamba-index:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-index sambamba/index.d

sambamba-sort:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-sort sambamba/sort.d

sambamba-view:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-view sambamba/view.d

sambamba-slice:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-slice sambamba/slice.d

sambamba-markdup:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -version=standalone -ofbuild/sambamba-markdup sambamba/markdup.d

.PHONY: clean

clean:
	rm -rf build/
