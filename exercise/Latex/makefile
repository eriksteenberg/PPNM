CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm
CC = gcc
default:ex.pdf ex.tex data.txt pyxplot.pdf

ex.pdf: ex.tex pyxplot.pdf 
	pdflatex $<
	pdflatex $<
#	bibtex $<
#	pdflatex $<

data.txt:main
	./$< > $@

main:main.o ex.o

ex.o: ex.c
	$(CC) $(CFLAGS) -c ex.c

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

pyxplot.pdf: pyxplot.gpi data.txt
	pyxplot $<

.PHONEY:clean
clean:
	$(RM) ex.png main main.o ex.o data.txt  *.pdf *.log *.aux 
