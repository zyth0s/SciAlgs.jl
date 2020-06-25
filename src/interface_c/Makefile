
CC=gcc

CFLAGS=-c -Wall -fPIC

SOURCES=calc_mean.c
OBJECTS=$(SOURCES:.c=.o)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

lib: $(OBJECTS)
	$(CC) -shared -fPIC -o libmean.so $(OBJECTS)

clean:
	rm *.o *.so
