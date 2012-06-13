all: alignment-test

sam_alignment.d: sam_alignment.rl
	ragel sam_alignment.rl -D -G2

alignment-test: sam_alignment.d
	rdmd -I../ -unittest --main sam_alignment.d

clean:
	rm sam_alignment.d
