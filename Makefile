FILES=bamfile.d rangetransformer.d chunkinputstream.d bgzfrange.d \
	  utils/inputrangechunks.d samheader.d \
	  bindings.d

all: scaffolds
	dmd $(FILES) -oflibbam.so -O -release -inline -shared

debug:
	dmd $(FILES) -oflibbam.so -debug -g -unittest -shared

scaffolds:
	dmd samheader.d generate_scaffolds.d -ofgenerate_scaffolds -J.
	./generate_scaffolds

clean:
	rm *.o
