CFLAGS = -Wall -std=gnu11
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

result.txt: main
	./$< > $@
	cat result.txt

clean:
	$(RM) main.o main result.txt
