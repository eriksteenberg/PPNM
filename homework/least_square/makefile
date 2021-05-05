CC = gcc
CFLAGS = -Wall -O1 -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default:result.txt leastsquarefit.png

leastsquarefit.png: result.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "least square fit";\
	plot \
	"$<" index 0 using 1:2:3 with yerrorbars title "data points"\
	,"$<" index 1 using 1:2 with lines title "fit"\
	,"$<" index 2 using 1:2 with lines title "fit + c0 and c1 error"\
	,"$<" index 3 using 1:2 with lines title "fit - c0 and c1 error"\
		' |gnuplot

result.txt: main
	./$< > $@
	cat result.txt

main: main.o functions.o functions.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@

functions.o: functions.c
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	$(RM) main.o main functions.o result.txt *.png 




