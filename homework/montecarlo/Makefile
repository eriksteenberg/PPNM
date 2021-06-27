CFLAGS = -Wall -O -std=gnu11
#CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += -lm#$(shell /usr/bin/gsl-config --libs)
CC = gcc

default: result.txt plot.png
	
result.txt: main
	./$< > $@
	cat result.txt

plot.png: plot.data.txt 
	echo '\
	set terminal png;\
	set log y ;\
	set xlabel "number of points N";\
	set ylabel "Actual error";\
	set title "Comparison of quasi-random and pseudo-random monte Carlo integration";\
	plot \
	"$<" using 1:2 with linespoints title "quasi-random"\
	,"$<" using 1:3 with linespoints title "pseudo-random"\
		' | pyxplot

clean:
	$(RM) main *.txt *.gpi *.png

