CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


erf.pyxploy.png: data.txt tabulated.data.txt 
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" using 1:2 with line title "math erf(x)"\
	,"$<" using 1:3 title "gsl erf(x)"\
	,"$<" using 1:4 title "my erf(x)"\
	,"tabulated.data.txt" using 1:2 title "tabulated data"\
	'| tee test.gpi | pyxplot

test.gpi:erf.pyxplot.png
	touch $@

data.txt: main
	./$< > $@

main:main.o erf.o

.PHONEY:clean
clean:
	$(RM) erf.o main.o main data.txt *.png *.gpi

test:
	echo $(CFLAGS)
	echo $(LDLIBS)
