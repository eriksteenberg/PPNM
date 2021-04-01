CFLAGS = -Wall -O1 -std=gnu11 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)


main: main.o functions.o functions.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

main.o: main.c functions.h
	$(CC) $(CFLAGS) -c $< -o $@

functions.o: functions.c
	$(CC) $(CFLAGS) -c $< -o $@
