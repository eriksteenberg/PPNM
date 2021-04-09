CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = clang

default: test.png main result.txt

result.txt: main #more if additional information is needed
	./$< > $@
	cat result.txt

test.png: result.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "test--dy/dt";\
	plot \
	"$<" index 0 using 1:2 with line title "1 og 2"\
	' | tee test.gpi | pyxplot


main: main.o #functions.o functions.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c# functions.h
	$(CC) $(CFLAGS) -c $< -o $@

#functions.o: functions.c
#	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(RM) main.o main result.txt test.png test.gpi 
