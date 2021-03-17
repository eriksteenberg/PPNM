GFLAGS = -Wall -Werror -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default:out.txt spline.png

spline.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key bottom;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "spline";\
	plot \
	"$<" using 1:2 with point title "My spline of my function"\
	,"$<" using 1:3 with line title "my function"\
	'| tee test.gpi |pyxplot

out.txt:main
	./main > out.txt
	cat out.txt

main: main.o                           
	$(CC) -o main main.o $(LDLIBS) 

main.o: main.c                   
	$(CC) $(CFLAGS) -c main.c


.PHONEY:clean
clean:
	$(RM) main.o main out.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)                             
