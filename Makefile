CC=g++
FLAGS=-Wall -fopenmp -O2


all: clean clustering.cpp
	$(CC) clustering.cpp -o clustering $(FLAGS)

clean:
	rm -f *.out *.o *.swp ranking.cluster* *.matrix
