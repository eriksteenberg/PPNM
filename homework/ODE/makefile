CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

out.txt: main
	./$< > $@
	cat out.txt

main: main.o #functions.o functions.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c# functions.h
	$(CC) $(CFLAGS) -c $< -o $@

#functions.o: functions.c
#	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(RM) main.o main out.txt
