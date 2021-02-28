CFLAGS = -Wall -Werror 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS =  $(shell gsl-config --libs)

result.txt:main
	./$< > $@

main:main.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o main main.c $(LDLIBS)

.PHONEY:clean
clean:
	$(RM) result.txt main main.o
