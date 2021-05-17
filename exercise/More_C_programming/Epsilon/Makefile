CFLAGS = -std=gnu11 -O -fwrapv
LDLIBS = -lm
CC = gcc

out.txt:main
	./$< > out.txt
	cat out.txt

main: main.o equal.o name_digit.o
	$(CC) $(CFLAGS) main.o equal.o name_digit.o -o main

name_digit.o: name_digit.c
	$(CC) $(CFLAGS) -c $< -o $@
	
equal.o: equal.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY:clean
clean:
	$(RM) *.o main out.txt
