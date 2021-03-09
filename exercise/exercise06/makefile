CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default: out.txt my_erf.png bessel.txt

my_erf.png: erf.data.txt tabulated.erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "error function";\
	plot \
	"$<" using 1:2 with line title "my\_erf(x)"\
	,"tabulated.erf.txt" using 1:2 with point title "tabulated erf points"\
	'| tee test.gpi | pyxplot
out.txt: main
	./$< > $@
	cat out.txt

erf.data.txt: erf
	./$< > $@

bessel.txt: bessel
	./$< > $@

main: main.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o main main.c $(LDLIBS)

erf: erf.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o erf erf.c $(LDLIBS)

bessel: bessel.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o bessel bessel.c $(LDLIBS)

.PHONEY:clean
clean:
	$(RM) out* main erf *.png *.gpi bessel bessel.txt


.PHONEY:test
test:
	echo CFLAGS
	echo $(CFLAGS)
	echo $(LDLIBS)
