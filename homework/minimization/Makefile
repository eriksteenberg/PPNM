CC = gcc
CFLAGS =  -O1 -std=gnu11 #-Wall
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: result.txt Higgs.png plot.txt

Higgs.png: data.txt plot.txt #data.txt #Makefile
	echo '\
	set terminal png;\
	set key top right ;\
	set tics out;\
	set xlabel "Energy (GeV)";\
	set ylabel "cross section";\
	set title "Mass of the Higgs boson estimate";\
	plot \
	"data.txt" using 1:2:3 with yerrorbars title "DATA" \
	,"plot.txt" index 0 using 1:2 with lines title "fit with mass = 125.97218" \
	' | pyxplot

#	"plot.txt" index 0 using 1:2 with lines title "fit" \
	,"data.txt" using 1:2 with lines title "Data"\ 

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
	$(RM) *.o main result.txt plot.txt *.png
