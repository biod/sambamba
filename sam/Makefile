all: alignment-test

fastrecordparser: sam_alignment.rl
	ragel sam_alignment.rl -D -G2
	echo 'module sam.recordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp recordparser.d

recordparser.d: sam_alignment.rl
	ragel sam_alignment.rl -D
	./fix_static_const.sh
	echo 'module sam.recordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp recordparser.d

alignment-test: recordparser.d
	rdmd -I../ -unittest --main recordparser.d

clean:
	rm -f recordparser.d
