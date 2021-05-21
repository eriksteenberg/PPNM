CC = gcc
FLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

all: out_read_stdin.txt out.txt out_read.txt

out.txt: cmdline
	./cmdline 1 2 3 4 5 > out.txt

out_read_stdin.txt: read-stdin
	echo 1 2 3 4 5 | ./read-stdin > out_read_stdin.txt


out_read.txt: read-file input.txt
	./read-file input.txt out_read.txt

input.txt:
	seq 0 0.5 10 > input.txt

cmdline: cmdline.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)	

read-stdin: read-stdin.o	
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)	

read-file: read-file.o
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

read-file.o: read-file.c
	$(CC) $(CFLAGS) -c read-file.c -o read-file.o

cmdline.o: cmdline.c
	$(CC) $(CFLAGS) -c $< -o $@
	
read-stdin.o: read-stdin.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean: 
	$(RM) *.txt *.o cmdline read-stdin read-file
