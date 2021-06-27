CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: out.txt time.txt test.txt timeplot.png

test.txt: test
	./$< 2 >$@
	./$< 21 >>$@
	./$< 12 >>$@

time.txt: time
	./$< 5 >$@
	./$< 10 >$@
	./$< 15 >>$@
	./$< 20 >$@
	./$< 25 >>$@
	./$< 30 >>$@
	./$< 35 >>$@
	./$< 40 >>$@
	./$< 45 >>$@
	./$< 50 >>$@
	./$< 55 >>$@
	./$< 60 >>$@
	./$< 65 >>$@
	./$< 70 >>$@
	./$< 75 >>$@
	./$< 80 >>$@
	./$< 85 >>$@
	./$< 90 >>$@
	./$< 95 >>$@
	./$< 100 >>$@
	./$< 105 >>$@

out.txt: main
	./$< > $@
	cat out.txt


	#f(x) = c1*(x**2)
	#fit f() "$<" using 1:2 via c1
	#f(x) with blue lines title "Fit to c_n*x**2"\


timeplot.png: time.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "N - dimension of matrices";\
	set ylabel "time";\
	set title "Time usage of n-sized matrices";\
	f(x) = c1*(x**2);\
	fit f() withouterrors "$<" using 1:2 via c1;\
	plot \
	"$<" using 1:2 with points title "data"\
	,f(x) with color blue lines title r"fit of $$cN^2$$ "\
	' | tee test.gpi | pyxplot


test: test.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	
test.o: test.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@
	
time: time.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	
time.o: time.c functions.o functions.h
	$(CC) $(CFLAGS) -c $< -o $@

main: main.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(RM) *.o *.txt time test main out.txt  *.png *.gpi
