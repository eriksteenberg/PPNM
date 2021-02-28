CFLAGS = `gsl-config --cflags`
CFLAGS = $(shell gsl-config --cflags) -O -Wall -std=gnu11
LDLIBS = $(shell gsl-config --libs)

start: ln.gamma.pyxploy.png gamma.pyxploy.png erf.pyxploy.png main gamma.txt ln.gamma.txt erf.txt

ln.gamma.pyxploy.png: ln.gamma.txt tabulated.lngamma.txt
	echo '\
        set terminal png;\
        set output "$@";\
        set key top left;\
        set tics out;\
        set xlabel "x";\
        set ylabel "y";\
	set title "logarithm of gamma";\
	plot \
        "$<" using 1:2 with line title "math lgamma(x)"\
        ,"$<" using 1:3 title "gsl lgamma(x)"\
        ,"$<" using 1:4 title "my lgamma(x)"\
	,"tabulated.lngamma.txt" using 1:2 title " tabulated data"\
        '| tee test.gpi | pyxplot i

gamma.pyxploy.png: gamma.txt tabulated.gamma.txt
	echo '\
        set terminal png;\
	set output "$@";\
        set key top left;\
        set tics out;\
        set xlabel "x";\
        set ylabel "y";\
        set title "gamma function";\
        plot \
        "$<" using 1:2 with line title "math gamma(x)"\
        ,"$<" using 1:3 title "gsl gamma(x)"\
        ,"$<" using 1:4 title "my gamma(x)"\
	,"tabulated.gamma.txt" using 1:2 title " tabulated data"\
        '| tee test.gpi | pyxplot 


erf.pyxploy.png: erf.txt tabulated.erf.txt
	echo '\
        set terminal png;\
        set output "$@";\
        set key top left;\
        set tics out;\
        set xlabel "x";\
	set ylabel "y";\
        set title "error function";\
        plot \
        "$<" using 1:2 with line title "math erf(x)"\
	,"$<" using 1:3 title "gsl erf(x)"\
	,"$<" using 1:4 title "my erf(x)"\
	,"tabulated.erf.txt" using 1:2 title " tabulated data"\
       '| tee test.gpi | pyxplot 

gamma.txt ln.gamma.txt erf.txt: main
	./$< > $@ 

main: main.o function.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LDLIBS)


function.o: function.c # erf.c gamma.c #ln.gamma.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(LDLIBS)

.PHONEY:clean
clean:
	$(RM) *.o main erf.txt gamma.txt ln.gamma.txt *.png *.gpi

test:
	echo $(CFLAGS)
	echo $(LDLIBS)
