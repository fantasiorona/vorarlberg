/* ==============================================================================================
Multimodal Problem

Apply a parallel approach (mpich) to the multimodal problem of the Ackley�s function.
==============================================================================================*/
// TODO: Run more than one population in parallel and switch out genoms at certain points
// TODO: Run on cluster with mpich
// TODO: think about better crossover maybe?
// TODO: Mabye try different versions of mutate and init functions for different populations
//#define VERBOSE
#include "Population.h"

#include <chrono>
#define _USE_MATH_DEFINES
#include <math.h>
//#include <mpi.h>

std::function<void(std::vector<double> &)>
get_initialization_function(const std::vector<double> &lower_gene_bounds,
                            const std::vector<double> &upper_gene_bounds);
double evaluate_genotype(std::vector<double> &genes);
void mutate_genotype(std::vector<double> &genes, double mutationProbability,
                     Population<double> &population);
void crossover_genotypes(std::vector<double> &genes1, std::vector<double> &genes2,
                         Population<double> &population);

const int maxGenerations = 1000000;

// TODO: Vary this for different populations
const int populationSize = 20;

// TODO: Vary this for different populations
const float crossoverChance = 0.8;

// TODO: Vary this for different populations
const float mutationChance = 0.5;

// in what area we try to find the global minima
// assumed to always be a square (or higher dimensional equivalent) with center of (0/0)
const double bound_extents = 5.0;
// in how many dimensions we want to build the ackley function
const unsigned int ackley_parameter = 4;

// TODO: Vary this for different populations
// how much we want to mutate values by (adds a value between -max_mutation_step and
// +max_mutation_step to current value)
const double max_mutation_step = 1.0;

// MPI variables
const int populationAmnt = 4; // should be a multiple of the Amount of processes we use
int procId, procAmnt, mpiError;
// MPI_Status status;

int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> lower_gene_bounds(ackley_parameter, -bound_extents);
    std::vector<double> upper_gene_bounds(ackley_parameter, bound_extents);

    std::vector<Population<double> *> populations;
    for (int i = 0; i < populationAmnt; i++) {
        auto population =
            new Population<double>(lower_gene_bounds, upper_gene_bounds, maxGenerations,
                                   populationSize, crossoverChance, mutationChance);
        populations.push_back(population);
    }

    bool allFinished = false;

    while (!allFinished) {
        // Evolve all populations
        for (int i = 0; i < populationAmnt; i++) {
            populations[i]->evolveParallel(
                get_initialization_function(lower_gene_bounds, upper_gene_bounds),
                evaluate_genotype, mutate_genotype, crossover_genotypes, 0.0);
        }

        allFinished = true;
        for (int i = 0; i < populationAmnt; i++) {
            if (!populations[i]->isFinished) {
                allFinished = false;
            }

            // Crossover between the populations
            for (int j = 0; j < populationAmnt; j++) {
                if (i != j) {
                    populations[i]->ReplaceWorstGenotypes(populations[j]->getGenotypes());
                }
            }
        }
    }

    // Print the result
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count() / 1000000;
    int milliseconds = (duration.count() - (seconds * 1000000)) / 1000;
    int microseconds = (duration.count() - (seconds * 1000000) - (milliseconds * 1000));

    std::cout << "Settings: " << std::endl;
    std::cout << "     Max generations:    " << maxGenerations << std::endl;
    std::cout << "     Populations size:   " << populationSize << std::endl;
    std::cout << "     Crossover chance:   " << crossoverChance << std::endl;
    std::cout << "     Mutation chance:    " << mutationChance << std::endl;
    std::cout << "     Lower Gene Bounds: ";
    for (auto gene : lower_gene_bounds)
        std::cout << gene << " ";
    std::cout << std::endl;
    std::cout << "     Upper Gene Bounds: ";
    for (auto gene : upper_gene_bounds)
        std::cout << gene << " ";
    std::cout << std::endl;
    std::cout << "Took: " << seconds << "s " << milliseconds << "ms " << microseconds << "us"
              << std::endl
              << std::endl;

    for (int i = 0; i < populationAmnt; i++) {
        std::cout << "--------- POPULATION #" << i << " ---------" << std::endl << std::endl;
        Genotype<double> bestMember = populations[i]->GetBestGenotype();
        populations[i]->print_result();
    }
}

// TODO: Is this even how we should do this????
/*int MPImain(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> lower_gene_bounds(ackley_parameter, -bound_extents);
    std::vector<double> upper_gene_bounds(ackley_parameter, bound_extents);

    std::vector<Population<double>> populations;
    for (int i = 0; i < populationAmnt; i++) {
        Population<double> population(lower_gene_bounds, upper_gene_bounds, maxGenerations,
                                      populationSize, crossoverChance, mutationChance);
        populations.push_back(population);
    }

    // Initialize the MPI environment
    mpiError = MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procAmnt);
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    // Evolve the population
    for (int i = 0; i < populationAmnt / procAmnt; i++) {
        int index = procId + (i * procAmnt);
        populations[index].evolve(get_initialization_function(lower_gene_bounds, upper_gene_bounds),
                                  evaluate_genotype, mutate_genotype, crossover_genotypes, 0.0);
    }

    // at this point the evolutions should be calculated
    if (procId == 0) {
        for (int source = 1; source < populationAmnt - populationAmnt / procAmnt; source++) {
            double genesNum;
            mpiError = MPI_Recv(&genesNum, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,
&status);

            double* genes;
            mpiError = MPI_Recv(genes, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            std::vector<double> genesVector(genes, genes + genesNum);
            double fitness = evaluate_genotype(genesVector);
            // TODO: evaluate genotypes and switch them??

#ifdef VERBOSE
            std::cout << "Received " << val << std::endl;
#endif
        }
    } else {
        for (int i = 0; i < populationAmnt / procAmnt; i++) {
            int index = procId + (i * procAmnt);

            Genotype<double> bestMember = populations[index].GetBestGenotype();
            std::vector<double> genes = bestMember.genes;

            size_t size = genes.size();
            MPI_Send(&size, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&genes.data(), genes.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    int bestIndex = procId;
    if (procId == 0) {
        // Evolve the population
        for (int i = 0; i < populationAmnt / procAmnt; i++) {
            int index = procId + (i * procAmnt);
            // TODO: find best genotype handled by process 0 and set bestIndex to its index
            Genotype<double> bestMember = populations[index].GetBestGenotype();
        }
    }
    mpiError = MPI_Finalize();

    // Print the result
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count() / 1000000;
    int milliseconds = (duration.count() - (seconds * 1000000)) / 1000;
    int microseconds = (duration.count() - (seconds * 1000000) - (milliseconds * 1000));

    std::cout << "Settings: " << std::endl;
    std::cout << "     Max generations:    " << maxGenerations << std::endl;
    std::cout << "     Populations size:   " << populationSize << std::endl;
    std::cout << "     Crossover chance:   " << crossoverChance << std::endl;
    std::cout << "     Mutation chance:    " << mutationChance << std::endl;
    std::cout << "     Lower Gene Bounds: ";
    for (auto gene : lower_gene_bounds)
        std::cout << gene << " ";
    std::cout << std::endl;
    std::cout << "     Upper Gene Bounds: ";
    for (auto gene : upper_gene_bounds)
        std::cout << gene << " ";
    std::cout << std::endl;
    std::cout << "Took: " << seconds << "s " << milliseconds << "ms " << microseconds << "us"
              << std::endl
              << std::endl;
    Genotype<double> bestMember = populations[bestIndex].GetBestGenotype();
    populations[bestIndex].print_result();
}*/

std::function<void(std::vector<double> &)>
get_initialization_function(const std::vector<double> &lower_gene_bounds,
                            const std::vector<double> &upper_gene_bounds) {
    return [lower_gene_bounds, upper_gene_bounds](std::vector<double> &genes) {
        int variableAmount = genes.size();
        // Set the genes to random numbers
        for (int gene = 0; gene < variableAmount; gene++) {
            float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            // finding a random value between both bounds
            double value =
                lower_gene_bounds[gene] + (upper_gene_bounds[gene] - lower_gene_bounds[gene]) * r;
            genes[gene] = value;
#ifdef VERBOSE
            std::cout << genes[gene] << ", ";
#endif
        }
#ifdef VERBOSE
        std::cout << std::endl;
#endif
    };
}

double evaluate_genotype(std::vector<double> &genes) {
    double a = 20;
    double b = 0.2;
    double c = 2.0 * M_PI;
    double d = genes.size();

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int i = 0; i < d; i++) {
        sum1 += genes[i] * genes[i];
        sum2 += cos(c * (double)genes[i]);
    }

    // ackley function
    double result = -a * exp(-b * sqrt(sum1 / d)) - exp(sum2 / d) + a + exp(1.0);

    // returning minus absolute value because we know best value would be 0, so everything that is
    // bigger or smaller is equaly bad
    // minus because library assumes higer value is better fitness
    return -abs(result);
}

// swap genes randomly without any constraints
void mutate_genotype(std::vector<double> &genes, double mutationProbability,
                     Population<double> &population) {
    int variableAmount = genes.size();
    // iterating over each gene
    for (int j = 0; j < variableAmount; ++j) {
        // random chance to mutate gene
        double x = population.GetRandomNormalizedDouble();
        if (x < mutationProbability) {
            // adding a random value in mutation_step range to current value
            double value =
                genes[j] + (population.GetRandomNormalizedDouble() * 2.0 - 1.0) * max_mutation_step;
#ifdef VERBOSE
            std::cout << "Mutating gene of value " << genes[j] << " to value " << value
                      << std::endl;
#endif
            genes[j] = value;
        }
    }
}

void crossover_genotypes(std::vector<double> &genes1, std::vector<double> &genes2,
                         Population<double> &population) {

    size_t geneAmount = genes1.size();
    // crossover is a simple swap operation with a randomly defined cutoff point
    int cutoff_point = population.get_random_gene_index();
#ifdef VERBOSE
    std::cout << "Genes before Crossover: " << std::endl;
    std::cout << "    Gene1: ";
    for (auto gene : genes1)
        std::cout << gene << ", ";
    std::cout << std::endl;
    std::cout << "    Gene2: ";
    for (auto gene : genes2)
        std::cout << gene << ", ";
    std::cout << std::endl;
    std::cout << "Crossing over at point: " << cutoff_point << ": " << std::endl;
#endif
    for (int i = cutoff_point; i < geneAmount; ++i) {
        std::swap(genes1[i], genes2[i]);
    }

#ifdef VERBOSE
    std::cout << "Genes after Crossover: " << std::endl;
    std::cout << "    Gene1: ";
    for (auto gene : genes1)
        std::cout << gene << ", ";
    std::cout << std::endl;
    std::cout << "    Gene2: ";
    for (auto gene : genes2)
        std::cout << gene << ", ";
    std::cout << std::endl;
#endif
}
