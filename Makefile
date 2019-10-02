# the compiler: gcc for C program, define as g++ for C++
CC = gcc
H5CC = h5cc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS = -O2 -Wall

DEBUG_FLAGS = -g -O0

# library flags:
#  -lz links zlib
#  -lm links math library
LIBFLAGS = -lz -lm

# the build target executable:
wgsim: wgsim.c tree.h kseq.h xrand.h utils.h
	$(CC) $(CFLAGS) -o wgsim wgsim.c $(LIBFLAGS)

debug: wgsim.c tree.h kseq.h xrand.h utils.h
	$(CC) $(DEBUG_FLAGS) -o wgsim wgsim.c $(LIBFLAGS)

tree2dmat: tree2dmat.c tree.h utils.h
	$(H5CC) $(CFLAGS) -o tree2dmat tree2dmat.c $(LIBFLAGS)
