default: out.txt
	cat out.txt
out.txt: main
	./main > out.txt
main: main.o
main.o: main.c
clean: 
	$(RM) *.o out* main
