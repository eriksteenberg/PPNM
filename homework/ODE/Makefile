CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: test.png SIR.png

result.txt: main
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
	set title "uÂ´Â´=-u u(0)=0 uÂ´(0)=1";\
	plot \
	"$<" index 0 using 1:2 with line title "u"\
	,"$<" index 0 using 1:3 with line title "dudt"\
	' | tee test.gpi | pyxplot

SIR.png: result.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "time [days]";\
	set ylabel "Peple";\
	set title "The SIR-model of epidemic developement, loosely based on Covid-19 info and the danish population with one infected";\
	plot \
	"$<" index 1 using 1:2 with line title "Susceptible"\
	,"$<" index 1 using 1:3 with line title "Infectious"\
	,"$<" index 1 using 1:4 with line title "Removed"\
	' | tee sir.gpi | pyxplot

main: main.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	$(RM) main.o main result.txt *.png *.gpi 
