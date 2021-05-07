CC = gcc
CFLAGS =  -O1 -std=gnu11 #-Wall
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

result.txt: main
	./$< > $@
	cat result.txt

main: main.o functions.o functions.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@

functions.o: functions.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o main Results.txt
