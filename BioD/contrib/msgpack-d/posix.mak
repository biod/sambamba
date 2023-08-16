# build mode: 32bit or 64bit
MODEL ?= $(shell getconf LONG_BIT)

ifeq (,$(DMD))
	DMD := dmd
endif

LIB    = libmsgpack-d.a
DFLAGS = -Isrc -m$(MODEL) -d -w -dip25 -dip1000

ifeq (true, $(EnableReal))
	DFLAGS += -version=EnableReal
endif

ifeq ($(BUILD),debug)
	DFLAGS += -g -debug
else
	DFLAGS += -O -release -nofloat -inline -noboundscheck
endif

NAMES = attribute common package register unpacker buffer exception packer streaming_unpacker  value
FILES = $(addsuffix .d, $(NAMES))
SRCS  = $(addprefix src/msgpack/, $(FILES))

# DDoc
DOCS      = $(addsuffix .html, $(NAMES))
DOCDIR    = html
CANDYDOC  = $(addprefix html/candydoc/, candy.ddoc modules.ddoc)
DDOCFLAGS = -Dd$(DOCDIR) -c -o- -Isrc $(CANDYDOC)

target: doc $(LIB)

$(LIB):
	$(DMD) $(DFLAGS) -lib -of$(LIB) $(SRCS)

doc:
	$(DMD) $(DDOCFLAGS) $(SRCS)

clean:
	rm $(addprefix $(DOCDIR)/, $(DOCS)) $(LIB)

MAIN_FILE = "empty_msgpack_unittest.d"

unittest:
	echo 'import msgpack; void main(){}' > $(MAIN_FILE)
	$(DMD) $(DFLAGS) -unittest -of$(LIB) $(SRCS) -run $(MAIN_FILE)
	rm $(MAIN_FILE)

run_examples:
	echo example/* | xargs -n 1 $(DMD) src/msgpack/*.d $(DFLAGS) -Isrc -run
