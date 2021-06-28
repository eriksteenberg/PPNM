CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: time.txt  timeplot.png

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
	./$< 95 >>$@
	./$< 100 >>$@
	./$< 105 >>$@
	./$< 110 >>$@
	./$< 115 >>$@
	./$< 135 >>$@
	./$< 150 >>$@
	./$< 175 >>$@
	./$< 200 >>$@
	./$< 215 >>$@
	./$< 230 >>$@
	./$< 250 >>$@
	./$< 275 >>$@
	./$< 300 >>$@

# Best fit parameters were:
# # -------------------------
#
# c0 = -0.60901188
# c1 = 0.0065626538
# c2 = 6.6195201e-05


timeplot.png: time.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel r"$$N-$$dimension of matrices";\
	set ylabel r"Time$$[s]$$";\
	set title r"Time usage of $$N-$$sized matrices eigenvalue calculation using rank-2 updates with jacobi diagonalization";\
	f(x) = c0 + c1*x + c2*(x**2);\
	fit f() withouterrors "$<" using 1:2 via c0,c1,c2;\
	plot \
	"$<" using 1:2 with points title "data"\
	,f(x) with color blue lines title r"fit of $$c_0 + c_1*N + c_2*N^2$$ "\
	' | tee test.gpi | pyxplot

	
time: time.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)
	
time.o: time.c functions.o functions.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONEY: clean
clean:
	$(RM) *.o time.txt out.txt time test main out.txt  timeplot.png *.gpi
