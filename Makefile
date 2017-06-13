CC=g++
CFLAGS=-g -Wall -O2 -Wno-unused-function
HEADERS=kstring.h
OBJECTS=
LIBS=-lz -lm

all: computeNF

computeNF: computeNF.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $< -o $@ $(LIBS)
