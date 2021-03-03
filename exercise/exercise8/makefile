CC = gcc
CFLAGS = -Wall -Werror -O1 -std=gnu11
CFLAGS += -pthread -fopenmp
LDLIBS += -lm -lgomp
LDLIBS += -pthread
#CC = gcc

default: main out.txt

out.txt: main makefile
	./$< > $@
	cat out.txt

main: main.o
	$(CC) $(CFALGS) $< -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONEY:clean
.PHONEY:clean
clean:
	$(RM) main *.o *.txt

test:
	@echo $(CFLAGS) 
	@echo $(LDLIBS)
