CC = gcc
FLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

out.txt: cmdline read-stdin input.txt
	./read-stdin < input.txt > out.txt
input.txt:
	seq 0 0.5 10 > input.txt

cmdline: cmdline.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)	

read-stdin: read-stdin.o	
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)	

cmdline.o: cmdline.c
	$(CC) $(CFLAGS) -c $< -o $@
	
read-stdin.o: read-stdin.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean: 
	$(RM) out.txt cmdline.o cmdline input.txt read-stdin.o	read-stdin
