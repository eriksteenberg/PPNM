GS = -Wall -Werror -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


#out.txt:main
#	./main > out.txt
#	cat out.txt

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
