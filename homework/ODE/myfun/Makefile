CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: test.png 


test.png: result.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "gravity is cool";\
	plot \
	"$<" index 0 using 1:2 with line title "planet 1"\
	,"$<" index 0 using 3:4 with line title "planet 2"\
	,"$<" index 0 using 5:6 with line title "planet 3"\
	' | tee test.gpi | pyxplot
result.txt: main
	./$< > $@
	cat result.txt

main: main.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	$(RM) main.o main result.txt *.png *.gpi 
