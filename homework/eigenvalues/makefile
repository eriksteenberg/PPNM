CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: out.txt B_part.txt plot.png

plot.png: eigen.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "plot of the lowest eigenfunctions";\
	plot \
	"$<" index 0 using 1:2 with line title "0th order"\
	,"$<" index 1 using 1:2 with line title "1th order"\
	,"$<" index 2 using 1:2 with line title "2nd order"\
	' | tee test.gpi | pyxplot

out.txt: main
	./$< > $@
	cat out.txt

B_part.txt: hamilton
	./$< > $@
	cat B_part.txt

hamilton: hamilton.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

hamilton.o: hamilton.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@

functions.o: functions.c
	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(RM) *.o main out.txt B_part.txt hamilton *.png *.gpi
