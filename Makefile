D_COMPILER=dmd
D_FLAGS=-IBioD #-O -release -inline 

RDMD_FLAGS=--compiler=$(D_COMPILER) --build-only $(D_FLAGS)

all:
	mkdir -p build/
	rdmd $(RDMD_FLAGS) -ofbuild/sambamba main.d

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

.PHONY: clean

clean:
	rm -rf build/
