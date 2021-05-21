GFLAGS = -Wall -Werror -O1 -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default:out.txt spline.png quad.png

spline.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "spline of f(x) = 0.5*x";\
	plot \
	"$<" index 1 using 1:3 with line lw 8 title "linear interp integ"\
	,"$<" index 2 using 1:3 with line lw 5 title "gsl interp integrate"\
        ,"$<" index 0 using 1:2 with line lw 8 title "The function"\
        ,"$<" index 1 using 1:2 with line lw 5 title "linear interp"\
        ,"$<" index 2 using 1:2 with line lw 1 title "gsl interp"\
	' | tee test.gpi | pyxplot


quad.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "quad spline of f(x) = 0.5*x";\
	plot \
	"$<" index 3 using 1:2 with line lw 3 title "quad interp"\
	,"$<" index 3 using 1:3 with line lw 2 title "quad interp integ"\
	,"$<" index 4 using 1:2 with line lw 5 title "quad diff"\
	' | tee test.gpi | pyxplot

#"$<" index 3 using 1:3 with line lw 2 title "quad interp integ"\

out.txt:main
	./main > out.txt
	cat out.txt

main: main.o                           
	$(CC) -o main main.o $(LDLIBS) 

main.o: main.c                   
	$(CC) $(CFLAGS) -c main.c


.PHONEY:clean
clean:
	$(RM) main.o main out.txt test.gpi *.png

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)                             
