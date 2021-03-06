CFLAGS = -Wall -O -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)
CC = gcc

default: out.txt  time.txt  timeplot.png


#THis is how I calculated the run time of larger and large matrices eigenvalue calculations. 
#I think I would have done this another way had I been able to use nested functions, but I wasn't.
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


out.txt: main
	./$< > $@
	cat out.txt


#A fit is calculated using pyxplot
#the factors calculated had the following values:
# Best fit parameters were:
# # -------------------------
#
# c0 = 0.20295899
# c1 = -0.0069328925
# c2 = 5.68657212e-05

timeplot.png: time.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel r"$$N-$$dimension of matrices";\
	set ylabel r"Time$$[s]$$";\
	set title r"Time usage of $$N-$$sized matrices";\
	f(x) = c0 + c1*x + c2*(x**2);\
	fit f() withouterrors "$<" using 1:2 via c0,c1,c2;\
	plot \
	"$<" using 1:2 with points title "data"\
	,f(x) with color blue lines title r"fit of $$c_0+c_1*N + c_2*N^2$$ "\
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
	$(RM) *.o time.txt out.txt time test main out.txt  timeplot.png *.gpi
