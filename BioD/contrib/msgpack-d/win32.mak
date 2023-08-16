DMD     = dmd
LIB     = msgpack.lib
DFLAGS  = -O -release -inline -nofloat -w -d -Isrc
UDFLAGS = -w -g -debug -unittest

SRCS = src\msgpack.d

# DDoc
DOCDIR    = html
CANDYDOC  = html\candydoc\candy.ddoc html\candydoc\modules.ddoc
DDOCFLAGS = -Dd$(DOCDIR) -c -o- -Isrc $(CANDYDOC)

DOCS = $(DOCDIR)\msgpack.html

target: doc $(LIB)

$(LIB):
	$(DMD) $(DFLAGS) -lib -of$(LIB) $(SRCS)

doc:
	$(DMD) $(DDOCFLAGS) $(SRCS)

clean:
	rm $(DOCS) $(LIB)
