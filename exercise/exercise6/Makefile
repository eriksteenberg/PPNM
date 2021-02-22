CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

#erf.gnuplot.png: test.gpi
#	cat test.gpi | sed 's/pyxplot/gnuplot/' | gnuplot

erf.pyxplot.png: data.txt Makefile tabulated.data.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "homemade Erf"\
	,"tabulated.data.txt" using 1:2 with points title "tabulated data"\
	' | tee test.gpi | pyxplot

test.gpi:erf.pyxplot.png

data.txt: main
	./$< > $@

main: main.o erf.o

.PHONEY:clean
	clean: 
		$(RM) erf.o main main.o data.txt test *.png
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
