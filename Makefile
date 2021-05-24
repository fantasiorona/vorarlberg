CXX=mpicxx.mpich
CXXFLAGS=-O3 -Wall -std=c++11

run: mm.out
	mpiexec -np 4 ./mm.out

run-cluster: mm.out
	srun -n 64 --mpi=pmi2 ./mm.out

a.out:
	$(CXX) $(CXXFLAGS) src/mm.cpp src/Population.h src/Genotype.h
	
mm.out: a.out
	mv a.out mm.out
