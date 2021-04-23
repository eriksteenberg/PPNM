CFLAGS = -Wall -O -std=gnu11
#CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = -lm# $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: result.txt

result.txt: main
	./$< > $@
	cat result.txt
clean:
	$(RM) main *.txt
