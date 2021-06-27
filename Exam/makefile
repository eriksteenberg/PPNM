CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: out.txt

out.txt: test
	./$< 2 >$@
	./$< 21 >>$@
	./$< 12 >>$@

#	./$< > $@
#	cat out.txt

test: test.o functions.o 
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	
test.o: test.c 
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(RM) *.o test main out.txt  *.png *.gpi
