CC = gcc
CFLAGS =  -O1 -std=gnu11 -Wall
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: result.txt output.txt plot.png

plot.png: output.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set title "Neural networks interpolarization";\
	plot \
	"$<" index 0 using 1:2 with lines title "interpolated data"\
	,"$<" index 1 using 1:2 with points title "data"\
	' | pyxplot

output.txt: main
	cat result

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
	$(RM) *.txt *.o main Results.txt
