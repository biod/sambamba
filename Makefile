FILES=bamreader.d rangetransformer.d baminputstream.d bgzfrange.d

all:
	dmd $(FILES) -ofbamreader -O -release -inline

debug:
	dmd $(FILES) -ofbamreader -debug -g

ldc2:
	ldc2 $(FILES) -ofbamreader -O5 -release
