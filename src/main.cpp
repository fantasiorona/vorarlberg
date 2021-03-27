#include "Genotype.h"
#include "Population.h"
#include <vector>

int main() {
    Population<double> population("input.txt", 3, 1000, 50, 0.8, 0.15);

    population.evolve([](std::vector<double> genes) {
        return (genes[1] * genes[1]) - (genes[1] * genes[2]) + genes[3];
    });
}
