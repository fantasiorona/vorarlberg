CC=g++
CFLAGS=-O3 -Wall

test: tpp.out
	./tpp.out 0 1 2 4 6 9 12 15 18 19 -i 1000 -c

tpp.out:
	$(CC) src/tpp.cpp src/Population.h src/Genotype.h -o tpp.out
