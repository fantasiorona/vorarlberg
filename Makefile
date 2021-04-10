CC=g++
CFLAGS=-O3 -Wall

test: ms.out
	./ms.out 1 && sleep 2 && ./ms.out 2 && sleep 2 && ./ms.out 3

ms.out:
	$(CC) src/main.cpp src/Population.h src/Genotype.h -o ms.out
