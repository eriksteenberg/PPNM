CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)


erf.pyxploy.png: data.txt tabulated.data.txt 
	set terminal png;\
	set output "data.txt";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	"$<" using 1:2\
	,"$<"using 1:3\
	| tee test.gpi | pyxplot



test.gpi:erf.pyxplot.png

data.txt: main
	./$< > $@

main:main.o erf.o

#main.o: erf.c main.c

#.PHONEY:clean:
#clean:
#	$(RM) erf.o main.o main #data.txt test *.png

test:
	echo $(CFLAGS)
	echo $(LDLIBS)
