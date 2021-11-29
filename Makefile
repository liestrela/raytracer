CC=g++
CC_FLAGS=-fopenmp -Wall -Wextra -pedantic -std=c++20 -O3

all: main.o ray clean

ray: main.o
	$(CC) -o ray main.o -fopenmp -lgomp

main.o: main.cc
	$(CC) -c -o $@ $< $(CC_FLAGS) -isystem./include

debug: main.cc
	$(CC) -c -g -o main.o main.cc $(CC_FLAGS)
	make ray
	make clean

clean:
	rm -f main.o