#pragma once

#include "Genotype.h"
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <string>

template <class T>
class Population {
  public:
    std::vector<Genotype<T>> genotypes;
    std::vector<Genotype<T>> new_genotypes;
    Genotype<T> current_best_genotype;

    /// Upper and lower bounds for the gene variables
    std::vector<T> upper_gene_bound;
    std::vector<T> lower_gene_bound;

    unsigned int variable_count;
    unsigned int max_generations;
    unsigned int size;

    double crossover_probability;
    double mutation_probability;

    std::random_device random_device;
    std::mt19937 mersenne_twister_engine;
    std::uniform_real_distribution<double> random_normalized_double;
    std::uniform_int_distribution<int> random_gene_index;

    Population(std::string filepath, unsigned int variable_count, unsigned int max_generations,
               unsigned int size, double crossover_probability, double mutation_probability)
        : variable_count(variable_count), max_generations(max_generations), size(size),
          crossover_probability(crossover_probability), mutation_probability(mutation_probability),
          mersenne_twister_engine(random_device()), random_normalized_double(0.0, 1.0),
          random_gene_index(0, variable_count - 1) {
        genotypes.resize(size);
        new_genotypes.resize(size);
        lower_gene_bound.resize(variable_count);
        upper_gene_bound.resize(variable_count);

        parse_bounds_file(filepath);
        initialize_genotypes();
    }

    void parse_bounds_file(std::string filepath) {
        std::ifstream input;
        input.open(filepath);

        if (!input) {
            std::cerr << "Cannot open input file!" << std::endl;
            exit(1);
        }

        // Parse and set bounds
        for (int i = 0; i < variable_count; i++) {
            double lower_bound;
            double upper_bound;

            input >> lower_bound >> upper_bound;

            lower_gene_bound[i] = lower_bound;
            upper_gene_bound[i] = upper_bound;
        }

        input.close();
    }

    void initialize_genotypes() {
        for (Genotype<T> &genotype : genotypes) {
            genotype.set_variable_count(variable_count);

            // Set the genes to random numbers
            for (int gene = 0; gene < variable_count; gene++) {
                genotype.genes[gene] = get_random_gene_value(gene);
            }
        }
    }

    double get_random_normalized_double() {
        return random_normalized_double(mersenne_twister_engine);
    }

    int get_random_gene_index() {
        return random_gene_index(mersenne_twister_engine);
    }

    int get_random_gene_value(unsigned int index) {
        std::uniform_int_distribution<int> random_gene_value(lower_gene_bound[index],
                                                             upper_gene_bound[index]);

        return random_gene_value(mersenne_twister_engine);
    }

    /// Iterates through the population and crosses over randomly selected pairs.
    void crossover_population() {
        int current_match;
        int match_count = 0;

        for (int current_genotype = 0; current_genotype < size; ++current_genotype) {
            // Pick a random number for testing against the pre-defined crossover probability
            double x = get_random_normalized_double();

            if (x < crossover_probability) {
                // This genotype will be involved in a crossover
                ++match_count;

                if (match_count % 2 == 0) {
                    crossover(current_match, current_genotype);
                } else {
                    current_match = current_genotype;
                }
            }
        }
    }

    /// Performs a crossover on the two given parents.
    void crossover(int one, int two) {
        // Randomly select the point until which the crossover will be performed
        int cutoff_point = get_random_gene_index();

        for (int i = 0; i < cutoff_point; i++) {
            std::swap(genotypes[one].genes[i], genotypes[two].genes[i]);
        }
    }

    // Elitist function: find the best and worse genotypes, remember the best, and replace the worst
    // by the previous best if the current best is worse than the previous best.
    void elitist() {
        int index_of_best;
        int index_of_worst;

        double best_fitness = genotypes[0].fitness;
        double worst_fitness = genotypes[0].fitness;

        // Find the best and the worst genotypes
        for (int i = 0; i < size - 1; ++i) {
            if (genotypes[i + 1].fitness < genotypes[i].fitness) {
                if (genotypes[i].fitness >= best_fitness) {
                    best_fitness = genotypes[i].fitness;
                    index_of_best = i;
                }
                if (genotypes[i + 1].fitness <= worst_fitness) {
                    worst_fitness = genotypes[i + 1].fitness;
                    index_of_worst = i + 1;
                }
            } else {
                if (genotypes[i].fitness <= worst_fitness) {
                    worst_fitness = genotypes[i].fitness;
                    index_of_worst = i;
                }
                if (best_fitness <= genotypes[i + 1].fitness) {
                    best_fitness = genotypes[i + 1].fitness;
                    index_of_best = i + 1;
                }
            }
        }

        if (best_fitness >= current_best_genotype.fitness) {
            // If the best fitness here is better than the previously best fitness, remember this as
            // the current best
            current_best_genotype = genotypes[index_of_best];
        } else {
            // If we got worse in this iteration, replace the current worst genotype with the
            // previously best one
            genotypes[index_of_worst] = current_best_genotype;
        }
    }

    /// Evaluates the fitness of all genotypes using the given evaluation function, which takes a
    /// gene as a parameter and returns the fitness. The results are saved in the `fitness` field of
    /// each genotype.
    void evaluate_all_fitnesses(std::function<double(std::vector<T>)> evaluation_function) {
        for (Genotype<T> &genotype : genotypes) {
            genotype.fitness = evaluation_function(genotype.genes);
        }
    }

    void initialize(std::string filename, int &seed) {
        // TODO: Implement
    }

    /// Returns the genotype with the highest fitness.
    void remember_best_genotype() {
        for (Genotype<T> &genotype : genotypes) {
            if (genotype.fitness > current_best_genotype.fitness) {
                // TODO: We might want to avoid the overhead of copying all the time
                current_best_genotype = genotype;
            }
        }
    }

    // Mutate random genes of all genotypes.
    void mutate_population() {
        // Iterate over all genes in every genotype
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < variable_count; j++) {
                // Draw a random number; if it is below the mutation probability, replace the gene
                // with a random new one
                double x = get_random_normalized_double();

                if (x < mutation_probability) {
                    genotypes[i].genes[j] = get_random_gene_value(j);
                }
            }
        }
    }

    // Print some metrics to the console.
    void report(int generation) {
        if (generation == 0) {
            std::cout << "\n";
            std::cout << "  Generation       Best            Average       Standard \n";
            std::cout << "  number           value           fitness       deviation \n";
            std::cout << "\n";
        }

        double sum = 0.0;
        double sum_square = 0.0;

        for (int i = 0; i < size; i++) {
            sum = sum + genotypes[i].fitness;
            sum_square = sum_square + genotypes[i].fitness * genotypes[i].fitness;
        }

        double avg = sum / (double)size;
        double square_sum = avg * avg * size;
        double stddev = sqrt((sum_square - square_sum) / (size - 1));
        double best_val = current_best_genotype.fitness;

        std::cout << "  " << std::setw(8) << generation << "  " << std::setw(14) << best_val << "  "
                  << std::setw(14) << avg << "  " << std::setw(14) << stddev << "\n";
    }

    // Selector function which generates a new population by keeping random good members of the
    // current one
    void create_new_population() {
        // Find the total fitness of the population
        double fitness_sum = 0.0;
        for (int i = 0; i < size; i++) {
            fitness_sum = fitness_sum + genotypes[i].fitness;
        }

        // Calculate the relative fitness of each member and save it
        for (int i = 0; i < size; i++) {
            genotypes[i].relative_fitness = genotypes[i].fitness / fitness_sum;
        }

        // Calculate the cumulative fitness of each member and save it
        genotypes[0].cumulative_fitness = genotypes[0].relative_fitness;
        for (int i = 1; i < size; i++) {
            genotypes[i].cumulative_fitness =
                genotypes[i - 1].cumulative_fitness + genotypes[i].relative_fitness;
        }

        // Select survivors using the cumulative fitnesses
        for (int i = 0; i < size; i++) {
            double p = get_random_normalized_double();

            if (genotypes[0].cumulative_fitness > p) {
                new_genotypes[i] = genotypes[0];
            } else {
                for (int j = 0; j < size - 1; j++) {
                    // If this genotypes' cumulative fitness is small, but the next one's is large,
                    // that means the next one caused a significant increase --> add it to the new
                    // genotypes
                    if (genotypes[j].cumulative_fitness <= p &&
                        genotypes[j + 1].cumulative_fitness > p) {
                        new_genotypes[i] = genotypes[j + 1];
                    }
                }
            }
        }

        // Overwrite the old population with the new one
        for (int i = 0; i < size; i++) {
            genotypes[i] = new_genotypes[i];
        }
    }

    void evolve(std::function<double(std::vector<T>)> evaluation_function) {
        // Initial evaluation (this is called within the generation functions later)
        evaluate_all_fitnesses(evaluation_function);
        remember_best_genotype();

        for (int generation = 0; generation < max_generations; generation++) {
            create_new_population();
            crossover_population();
            mutate_population();
            report(generation);
            evaluate_all_fitnesses(evaluation_function);
            elitist();
        }

        // Report
        std::cout << "\n";
        std::cout << "  Best member after " << max_generations << " generations:\n";
        std::cout << "\n";

        for (int i = 0; i < variable_count; i++) {
            std::cout << "  var(" << i << ") = " << current_best_genotype.genes[i] << "\n";
        }

        std::cout << "\n";
        std::cout << "  Best fitness = " << current_best_genotype.fitness << "\n";

        // Done!
        std::cout << "\n";
        std::cout << "SIMPLE_GA:\n";
        std::cout << "  Normal end of execution.\n";
        std::cout << "\n";
    }
};