CC=/usr/bin/g++
#CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=bmean

all: $(EXEC)

bmean:   main.o  utils.o xor.o xxhash.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp   utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

xhash.o: xxhash.c xxhash.h
	$(CC) -o $@ -c $< $(CFLAGS)

xor.o: xor.cpp xor.h xxhash.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
