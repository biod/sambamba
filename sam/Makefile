all: alignment-test

sam_alignment.d: sam_alignment.rl
	ragel sam_alignment.rl -D -G2
	echo 'module sam.sam_alignment;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	mv .sam_alignment.d.tmp sam_alignment.d

alignment-test: sam_alignment.d
	rdmd -I../ -unittest --main sam_alignment.d

clean:
	rm sam_alignment.d
