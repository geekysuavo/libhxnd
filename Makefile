
# CC: compiler binary filename.
CC=gcc

# CFLAGS: compilation flags.
CFLAGS=-g -O2 -std=c99 -Wall -Wformat -I. -fopenmp
CFLAGS+= $(shell pkg-config --cflags gtk+-3.0)

# LIBS, GLIBS: linkage flags.
LIBS=-lm -lcblas -llapack
GLIBS=$(shell pkg-config --libs gtk+-3.0)

# LIBSRC: library source basenames: hypercomplex data structures.
LIBSRC=hx-algebra hx-scalar hx-index hx-array hx-array-mem hx-array-io
LIBSRC+= hx-array-rawio hx-array-topo hx-array-resize hx-array-slice
LIBSRC+= hx-array-tile hx-cmp hx-arith hx-blas hx-phasor hx-fourier
LIBSRC+= hx-window hx-ist

# LIBSRC: library source basenames: auxiliary library routines.
LIBSRC+= trace opts str bytes

# LIBSRC: library source basenames: core nmr datum structure.
LIBSRC+= nmr-datum nmr-datum-mem nmr-datum-io nmr-datum-array
LIBSRC+= nmr-datum-type nmr-datum-dims

# LIBSRC: library source basenames: nmr data formats.
LIBSRC+= nmr-hxnd nmr-text nmr-bruker nmr-varian nmr-pipe nmr-ucsf nmr-nv
LIBSRC+= nmr-rnmrtk

# LIBSRC: library source basenames: multivariate math functions.
LIBSRC+= mx-scaling

# LIBSRC: library source basenames: processing functions.
LIBSRC+= fn fn-args fn-list fn-abs fn-add fn-complex fn-cut fn-fft fn-ht fn-ist
LIBSRC+= fn-multiply fn-phase fn-real fn-resize fn-shift fn-window fn-zerofill

# GUISRC: graphical interface source basenames.
GUISRC=ghx-proc-main

# LIBOBJ, GUIOBJ: library and graphical interface object filenames.
LIBOBJ=$(addprefix libhxnd/,$(addsuffix .o,$(LIBSRC)))
GUIOBJ=$(addprefix ghx/,$(addsuffix .o,$(GUISRC)))

# BIN, BINBIN, BINOBJ: binary source, output and object filenames.
BIN=hx ghx
BINBIN=$(addprefix bin/,$(BIN))
BINOBJ=$(addprefix bin/,$(addsuffix .o,$(BIN)))

# OBJ: all object files that need compilation from source.
OBJ=$(LIBOBJ) $(GUIOBJ)

# registered suffixes for make rules.
.SUFFIXES: .c .o

# all: global, default compilation rule.
all: $(OBJ) $(BINBIN)

# bin/hx: binary linkage target for command-line hx multi-tool.
bin/hx: $(LIBOBJ) bin/hx.o
	@echo " LD $@"
	@$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

# bin/ghx: binary linkage target for graphical hx application.
bin/ghx: $(OBJ) bin/ghx.o
	@echo " LD $@"
	@$(CC) $(CFLAGS) $^ -o $@ $(LIBS) $(GLIBS)

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

# diff: target to diff all tracked changes in the git repo.
diff:
	@echo " GIT diff"
	@git diff --color

# git: target to commit and push all changes to github.
git: clean
	@echo " GIT commit"
	@git commit -a
	@echo " GIT push"
	@git push

