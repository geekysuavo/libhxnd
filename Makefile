
# CC, CFLAGS, LD: compilation and linkage flags.
CC=gcc
CFLAGS=-g -O2 -std=c99 -Wall -Wformat -I.
LIBS=-lm

# LIBSRC: library source basenames.
LIBSRC=hx-algebra hx-scalar hx-index hx-array hx-cmp hx-arith hx-fourier
LIBSRC+= trace opts str bytes nmr-datum nmr-bruker nmr-varian nmr-pipe
LIBSRC+= fn fn-fft fn-ht

# LIBOBJ: library object filenames.
LIBOBJ=$(addprefix libhxnd/,$(addsuffix .o,$(LIBSRC)))

# BIN, BINBIN, BINOBJ: binary source, output and object filenames.
BIN=hx
BINBIN=$(addprefix bin/,$(BIN))
BINOBJ=$(addprefix bin/,$(addsuffix .o,$(BIN)))

# OBJ: all object files that need compilation from source.
OBJ=$(LIBOBJ)

# registered suffixes for make rules.
.SUFFIXES: .c .o

# all: global, default compilation rule.
all: $(OBJ) $(BINBIN)

bin/hx: $(OBJ) bin/hx.o
	@echo " LD $@"
	@$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

clean:
	@echo " CLEAN"
	@rm -rf $(OBJ) $(BINOBJ) $(BINBIN)

again: clean all

fixme:
	@echo " FIXME"
	@grep -RHni fixme hxnd/*.h libhxnd/*.c bin/*.[ch]

linecount:
	@echo " WC"
	@wc -l hxnd/*.h libhxnd/*.c bin/*.[ch]

dist: clean
	@isodate=$$(date +%Y%m%d) && \
	 projdir=$$(pwd) && \
	 projdir=$$(basename $${projdir}) && \
	 echo " DIST $${isodate}" && \
	 cd .. && \
	 rm -rf $${projdir}-$${isodate}.tar.gz && \
	 tar cf $${projdir}-$${isodate}.tar $${projdir}/ && \
	 gzip -9 $${projdir}-$${isodate}.tar && \
	 cd $${projdir} >/dev/null

