CC=g++
CFLAGS=-O3 -Wall

test: ms.out
	./ms.out 1

ms.out:
	$(CC) src/main.cpp src/Population.h src/Genotype.h -o ms.out
