CFLAGS = -Wall -Werror 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS =  $(shell gsl-config --libs)
CC = gcc

default: result.txt erf.plot.png
	cat result.txt

erf.plot.png: erf.data.txt tabulated.erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "error function";\
	plot \
	"$<" using 1:2 with line title "my error function"\
	,"tabulated.erf.txt" using 1:2 with points title "tabulated data"\
	'| tee test.gpi | pyxplot

result.txt:main
	./$< > $@

erf.data.txt:main
	./$< > $@



main:main.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o main main.c $(LDLIBS)

.PHONEY:clean
clean:
	$(RM) result.txt main main.o erf.data.txt test.gpi *.png


#erf.data.txt  erf.plot.png  main.c  makefile  tabulated.erf.txt  test.gpi

