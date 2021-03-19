GFLAGS = -Wall -Werror -O1 -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


default:out.txt spline.png

spline.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "spline";\
	plot \
         "$<" index 1 using 1:3 with line ps 0.45 title "interpinteg"\
	,"$<" index 2 using 1:3 with line ps 0.35 title "gsl interp integrate"\
	,"$<" index 3 using 1:2 with line ps 0.15 title "quad interp"\
        ,"$<" index 0 using 1:2 with line ps 0.45 title "The function"\
        ,"$<" index 1 using 1:2 with line ps 0.35 title "my interp"\
        ,"$<" index 2 using 1:2 with line ps 0.35 title "gsl interp"\
        ,"$<" index 2 using 1:3 with line ps 0.35 title "gsl interp integrate"\
	' | tee test.gpi | pyxplot

out.txt:main
	./main > out.txt
	cat out.txt

main: main.o                           
	$(CC) -o main main.o $(LDLIBS) 

main.o: main.c                   
	$(CC) $(CFLAGS) -c main.c


.PHONEY:clean
clean:
	$(RM) main.o main out.txt test.gpi

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)                             
