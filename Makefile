
# CC, CFLAGS, LD: compilation and linkage flags.
CC=gcc
CFLAGS=-g -O2 -std=c99 -Wall -Wformat -I. -fopenmp
LIBS=-lm

# LIBSRC: library source basenames: hypercomplex data structures.
LIBSRC=hx-algebra hx-scalar hx-index hx-array hx-array-resize hx-cmp
LIBSRC+= hx-arith hx-phasor hx-fourier hx-window hx-ist

# LIBSRC: library source basenames: auxiliary library routines.
LIBSRC+= trace opts str bytes

# LIBSRC: library source basenames: nmr data formats.
LIBSRC+= nmr-datum nmr-bruker nmr-varian nmr-pipe nmr-ucsf

# LIBSRC: library source basenames: processing functions.
LIBSRC+= fn fn-abs fn-add fn-complex fn-cut fn-fft fn-ht fn-ist fn-phase
LIBSRC+= fn-real fn-resize fn-scale fn-shift fn-window fn-zerofill

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

# bin/hx: binary linkage target for command-line hx multi-tool.
bin/hx: $(OBJ) bin/hx.o
	@echo " LD $@"
	@$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

# .c.o: general compilation target for C source files.
.c.o:
	@echo " CC $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

# clean: remove all generated object code and binaries.
clean:
	@echo " CLEAN"
	@rm -rf $(OBJ) $(BINOBJ) $(BINBIN)

# again: quick full recompilation target.
again: clean all

# fixme: target to search all source files for 'fix me' statements.
fixme:
	@echo " FIXME"
	@grep \
	   --recursive --with-filename \
	   --line-number --ignore-case --color \
	 fixme hxnd/*.h libhxnd/*.c bin/*.[ch] man/*.[0-9] || \
	 echo " No statements found"

# linecount: target to count lines of all C source files and headers.
linecount:
	@echo " WC"
	@wc -l hxnd/*.h libhxnd/*.c bin/*.[ch]

# dist: target to generate a source tarball.
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

# git: target to commit and push all changes to github.
git: clean
	@echo " GIT commit"
	@git commit -a
	@echo " GIT push"
	@git push

