CFLAGS = -o -w -std=gnu99
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

all: output.txt 
	cat $<

output.txt: main
	./$< > $@ 

main: main.o
	$(CC) $(CFALGS) $< -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	$(RM) main output.txt main.o 
