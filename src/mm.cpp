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
#include <cstring>
#include <math.h>
#include <mpi.h>

std::function<void(std::vector<double> &)>
get_initialization_function(const std::vector<double> &lower_gene_bounds,
                            const std::vector<double> &upper_gene_bounds);
double evaluate_genotype(std::vector<double> &genes);
void mutate_genotype(std::vector<double> &genes, double mutationProbability,
                     Population<double> &population);
void crossover_genotypes(std::vector<double> &genes1, std::vector<double> &genes2,
                         Population<double> &population);

const int maxGenerations = 10000000;

// in what area we try to find the global minima
// assumed to always be a square (or higher dimensional equivalent) with center of (0/0)
const double bound_extents = 5.0;
// in how many dimensions we want to build the ackley function
const unsigned int parameterCount = 4;

// TODO: Vary this for different populations
// how much we want to mutate values by (adds a value between -max_mutation_step and
// +max_mutation_step to current value)
const double max_mutation_step = 1.0;

struct GenotypeData {
    double fitness;
    double relative_fitness;
    double cumulative_fitness;

    double genes[parameterCount];
};

int main(int argc, char **argv) {
    int procId, procAmnt, mpiError;
    MPI_Status status = {};

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &procAmnt);
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    std::cout << "\nCalculation started on process #" << procId << " ..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> lower_gene_bounds(parameterCount, -bound_extents);
    std::vector<double> upper_gene_bounds(parameterCount, bound_extents);

    int leftNbr = procId - 1;
    int rightNbr = procId + 1;
    if (leftNbr < 0) {
        leftNbr = procAmnt - 1;
    }
    if (rightNbr >= procAmnt) {
        rightNbr = 0;
    }

    // Parameters which are varied based on the process
    std::srand(procId);

    // I'd really like to make this parameter variable as well but I can't 🥲
    int populationSize = 30;

    // Range 0.5-0.9
    float crossoverChance = 0.5f + ((std::rand() % 40) / 100.f);

    // Range 0.5-0.9
    float mutationChance = 0.5f + ((std::rand() % 40) / 100.f);

#ifdef VERBOSE
    std::cout << "Proc ID: " << procId << " left: " << leftNbr << ", right: " << rightNbr
              << std::endl;
#endif

    auto population = new Population<double>(lower_gene_bounds, upper_gene_bounds, maxGenerations,
                                             populationSize, crossoverChance, mutationChance);

    auto data = new GenotypeData[populationSize];
    while (!population->isFinished) {
        // Evolve all populations
        population->evolveParallel(
            get_initialization_function(lower_gene_bounds, upper_gene_bounds), evaluate_genotype,
            mutate_genotype, crossover_genotypes, 0.0, procId);

        // Prepare data for sending
        for (int i = 0; i < population->getGenotypes().size(); i++) {
            auto &genotype = population->getGenotypes()[i];
            data->fitness = genotype.fitness;
            data->relative_fitness = genotype.relative_fitness;
            data->cumulative_fitness = genotype.cumulative_fitness;
            memcpy(data->genes, genotype.genes.data(), sizeof(double) * parameterCount);
            data->fitness = genotype.fitness;
        }

        MPI_Send(data, sizeof(GenotypeData) * populationSize, MPI_BYTE, rightNbr, 0,
                 MPI_COMM_WORLD);
        mpiError = MPI_Recv(data, sizeof(GenotypeData) * populationSize, MPI_BYTE, leftNbr, 0,
                            MPI_COMM_WORLD, &status);

        if (mpiError == MPI_SUCCESS) {
            // Create an std::vector<Genotype> from the received data
            int count;
            MPI_Get_count(&status, MPI_BYTE, &count);
            count /= sizeof(GenotypeData);

            std::vector<Genotype<double>> genotypes(count);

            for (int i = 0; i < genotypes.size(); i++) {
                genotypes[i].fitness = data->fitness;
                genotypes[i].relative_fitness = data->relative_fitness;
                genotypes[i].cumulative_fitness = data->cumulative_fitness;
                genotypes[i].genes.resize(parameterCount);
                memcpy(genotypes[i].genes.data(), data->genes, sizeof(double) * parameterCount);
            }
            population->ReplaceWorstGenotypes(genotypes);
        }
    }

    // Print the result
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    int seconds = duration.count() / 1000000;
    int milliseconds = (duration.count() - (seconds * 1000000)) / 1000;
    int microseconds = (duration.count() - (seconds * 1000000) - (milliseconds * 1000));

    Genotype<double> genotype();

    // Prepare data for sending
    double sendData, recvData, bestData;
    int bestIndex = 0;
    if (procId != 0) {
        std::vector<double> genes = population->GetBestGenotype().genes;
        sendData = evaluate_genotype(genes);

        MPI_Send(&sendData, 1, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD); // send data to process 0
    } else {
        std::vector<double> genes = population->GetBestGenotype().genes;
        bestData = evaluate_genotype(genes);

        for (int i = 1; i < procAmnt; i++) {
            mpiError = MPI_Recv(&recvData, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
            if (mpiError == MPI_SUCCESS) {
                if (recvData > bestData) {
                    bestData = recvData;
                    bestIndex = i;
                }
            }
        }
    }

    if (procId == 0) {
        std::cout << "\nSettings: " << std::endl;
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
    }

    population->print_result(procId);

    MPI_Bcast(&bestData, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestIndex, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::cout << "\nBest overall fitness: " << bestData << " from process #" << bestIndex
              << std::endl;

    mpiError = MPI_Finalize();
}

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

    // returning minus absolute value because we know best value would be 0, so everything that
    // is bigger or smaller is equaly bad minus because library assumes higer value is better
    // fitness
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
