CC=g++
FLAGS=-Wall


all: clean clustering.cpp
	$(CC) clustering.cpp -o clustering $(FLAGS)

clean:
	rm -f *.out *.o *.swp	
