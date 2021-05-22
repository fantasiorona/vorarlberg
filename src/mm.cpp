/* ==============================================================================================
Multimodal Problem

Apply a parallel approach (mpich) to the multimodal problem of the Ackley’s function.
==============================================================================================*/
#include "Population.h"

#include <chrono>

void initialze_genotype(std::vector<int> &genes);
double evaluate_genotype(std::vector<int> &genes);
void mutate_genotype(std::vector<int> &genes, double mutationProbability,
                     Population<int> &population);
void crossover_genotypes(std::vector<int> &genes1, std::vector<int> &genes2,
                         Population<int> &population);

// Declared here so they're re-used in each crossover sequence
std::vector<int> crossoverInverseSequence1;
std::vector<int> crossoverInverseSequence2;
std::vector<int> crossoverPositionSequence1;
std::vector<int> crossoverPositionSequence2;

// TODO: Adjust and create input file
const int maxGenerations = 1000000000;
const int populationSize = 20;
const float crossoverChance = 0.8;
const float mutationChance = 0.05;
const std::string file = "inputs/.....filename.....";

Population<int> population(file, maxGenerations, populationSize, crossoverChance, mutationChance);

int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();

    // adjusting size according to variable count
    int variableCount = population.GetVariableCount();
    crossoverInverseSequence1.resize(variableCount);
    crossoverInverseSequence2.resize(variableCount);
    crossoverPositionSequence1.resize(variableCount);
    crossoverPositionSequence2.resize(variableCount);

    // Evolve the population
    // TODO: implement evolve function without perfect fitness?
    population.evolve(initialze_genotype, evaluate_genotype, mutate_genotype, crossover_genotypes);

    // Print the result
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count() / 1000000;
    int milliseconds = (duration.count() - (seconds * 1000000)) / 1000;
    int microseconds = (duration.count() - (seconds * 1000000) - (milliseconds * 1000));

    std::cout << std::endl << std::endl << "Result for " << file << ":" << std::endl;
    std::cout << "Settings: " << std::endl;
    std::cout << "     Max generations:    " << maxGenerations << std::endl;
    std::cout << "     Populations size:   " << populationSize << std::endl;
    std::cout << "     Crossover chance:   " << crossoverChance << std::endl;
    std::cout << "     Mutation chance:    " << mutationChance << std::endl;
    std::cout << "Took: " << seconds << "s " << milliseconds << "ms " << microseconds << "us"
              << std::endl
              << std::endl;
    Genotype<int> bestMember = population.GetBestGenotype();
    population.print_result();
}

// TODO: check if this is correct
void initialze_genotype(std::vector<int> &genes) {
    int variableAmount = genes.size();
    std::vector<bool> usedValues(variableAmount, false);
    // Set the genes to random numbers
    for (int gene = 0; gene < variableAmount; gene++) {
        int value = population.GetRandomGeneValue(gene);
        // counting up if value is already in use, so we never get duplicate values
        while (usedValues[value - 1]) {
            ++value;
            if (value == variableAmount + 1) value = 1;
        }
        usedValues[value - 1] = true;
        genes[gene] = value;
    }
}

double evaluate_genotype(std::vector<int> &genes) {
    // TODO: implement ackley function
}

// swap genes randomly without any constraints
void mutate_genotype(std::vector<int> &genes, double mutationProbability,
                     Population<int> &population) {
    int variableAmount = genes.size();
    // iterating over each gene
    for (int j = 0; j < variableAmount; ++j) {
        // random chance to mutate gene
        double x = population.GetRandomNormalizedDouble();
        if (x < mutationProbability) {
            int swapIndex = population.GetRandomGeneValue(0) - 1;
            // make sure we dont swap in place
            while (swapIndex == j) {
                ++swapIndex;
                if (swapIndex == variableAmount + 1) swapIndex = 1;
            }
            std::swap(genes[j], genes[swapIndex]);
        }
    }
}

// TODO: check if this is correct
void crossover_genotypes(std::vector<int> &genes1, std::vector<int> &genes2,
                         Population<int> &population) {
    size_t geneAmount = genes1.size();

    // iterate over each value
    for (int i = 0; i < geneAmount; ++i) {
        // initialize to zero
        crossoverInverseSequence1[i] = 0;
        crossoverInverseSequence2[i] = 0;
        bool found1 = false;
        bool found2 = false;
        // for each value we start at the left
        for (int j = 0; j < geneAmount; ++j) {
            // and continue until we find the value
            if (!found1) {
                // if we pass a bigger value we increment
                if (genes1[j] > i + 1) {
                    crossoverInverseSequence1[i] = crossoverInverseSequence1[i] + 1;
                }
                // stop when we find the value
                if (genes1[j] == i + 1) {
                    found1 = true;
                }
            }
            // do same for the second parent
            if (!found2) {
                if (genes2[j] > i + 1) {
                    crossoverInverseSequence2[i] = crossoverInverseSequence2[i] + 1;
                }
                if (genes2[j] == i + 1) {
                    found2 = true;
                }
            }
        }
    }
    // crossover afterwards is a simple swap operation with a randomly defined cutoff point
    int cutoff_point = population.GetRandomGeneValue(0) - 1;
    for (int i = cutoff_point; i < geneAmount; ++i) {
        std::swap(crossoverInverseSequence1[i], crossoverInverseSequence2[i]);
    }
    // translate back to actual genes
    // we start with creating a position vector
    // this happens by counting greater or equal elements to the right of the current value
    // this is basically the inverse operation to the inverse sequence creation and gives the
    // information where each value is positioned in the gene vector
    for (int i = geneAmount - 1; i >= 0; --i) {
        // initialaze to current inverse value
        crossoverPositionSequence1[i] = crossoverInverseSequence1[i];
        crossoverPositionSequence2[i] = crossoverInverseSequence2[i];
        // for each value to the right
        for (int j = i + 1; j < geneAmount; ++j) {
            // increase if its bigger or equal
            if (crossoverPositionSequence1[j] >= crossoverPositionSequence1[i]) {
                crossoverPositionSequence1[j] = crossoverPositionSequence1[j] + 1;
            }
            // do the same for the other parent
            if (crossoverPositionSequence2[j] >= crossoverPositionSequence2[i]) {
                crossoverPositionSequence2[j] = crossoverPositionSequence2[j] + 1;
            }
        }
    }

    // actually create genes by using the position information to find the position of the value
    // add 1 because values are from 1-9 while vector indizes are from 0-8
    for (int i = 0; i < geneAmount; ++i) {
        genes1[crossoverPositionSequence1[i]] = i + 1;
        genes2[crossoverPositionSequence2[i]] = i + 1;
    }
}