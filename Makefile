# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS = -g -O2 -Wall

# library flags:
#  -lz links zlib
#  -lm links math library
LIBFLAGS = -lz -lm

# the build target executable:
wgsim: wgsim.c tree.h kseq.h
	$(CC) $(CFLAGS) -o wgsim wgsim.c $(LIBFLAGS)

