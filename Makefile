
all: wgsim.c kseq.h
	gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm 
