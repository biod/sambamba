all: fastrecordparser recordparser regionparser

.PHONY : fastrecordparser

.PHONY : recordparser

.PHONY : regionparser

fastrecordparser:
	ragel sam_alignment.rl -D -G2
	./workarounds/fix_switch_case_fallthrough.sh sam_alignment.d
	echo 'module bio.sam.utils.fastrecordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp fastrecordparser.d
	mv fastrecordparser.d ../bio/sam/utils/fastrecordparser.d

recordparser:
	ragel sam_alignment.rl -D
	./workarounds/fix_static_const.sh sam_alignment.d
	echo 'module bio.sam.utils.recordparser;' | cat - sam_alignment.d > .sam_alignment.d.tmp
	rm sam_alignment.d
	mv .sam_alignment.d.tmp recordparser.d
	mv recordparser.d ../bio/sam/utils/recordparser.d

regionparser:
	ragel region.rl -D
	./workarounds/fix_static_const.sh region.d
	mv region.d ../bio/core/region.d

clean:
	rm -f *parser.d region.d
