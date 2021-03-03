CFLAGS = -Wall -Werror
CFLAGS += $(shell gsl-config --cflags)
#CFLAGS += -pthread
LDLIBS += $(shell gsl-config --libs)
LDLIBS += -lm
CC = gcc

out.txt: main

main:main.o
	$(CC) $(CFALGS) $^ -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) main *.o *.txt

