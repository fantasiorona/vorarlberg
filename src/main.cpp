#include "Genotype.h"
#include "Population.h"
#include <vector>

int main() {
    // Create a new Population with genes consisting of three `double` variables
    Population<double> population("input.txt", 3, 1000, 50, 0.8, 0.15);

    // Evolve the population to maximize the given equation
    population.evolve([](std::vector<double> genes) {
        return (genes[0] * genes[0]) - (genes[0] * genes[1]) + genes[2];
    });

    // Print the result
    population.print_result();
}
