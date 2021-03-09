CC = gcc
CFLags = -std=gnu99 -o -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt
out.txt: math
	./math > out.txt
math: math.o
math.o:math.c

clean:
	$(RM) *.o out* math
