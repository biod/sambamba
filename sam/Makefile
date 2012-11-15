all: alignment-test

.PHONY : fastrecordparser

.PHONY : recordparser.d

fastrecordparser: 
	ragel sam_alignment.rl -D -G2
	../ragel_workarounds/fix_switch_case_fallthrough.sh sam_alignment.d
	echo 'module sam.recordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp recordparser.d

recordparser.d: 
	ragel sam_alignment.rl -D
	../ragel_workarounds/fix_static_const.sh sam_alignment.d
	echo 'module sam.recordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp recordparser.d

alignment-test: recordparser.d
	rdmd -I../ -unittest --main recordparser.d

clean:
	rm -f recordparser.d
