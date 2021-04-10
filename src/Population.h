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
    Population(std::string filepath, unsigned int max_generations, unsigned int size,
               double crossover_probability, double mutation_probability)
        : max_generations(max_generations), size(size),
          crossover_probability(crossover_probability), mutation_probability(mutation_probability),
          mersenne_twister_engine(random_device()), random_normalized_double(0.0, 1.0) {

        genotypes.resize(size);
        new_genotypes.resize(size);

        parse_bounds_file(filepath);
        // has to be initialized after parsing because variable_count only gets the correct value
        // after the file is read
        random_gene_index = std::uniform_int_distribution<int>(0, variable_count - 1);
    }

    // Evolve all generations and print intermediate results
    void evolve(std::function<void(std::vector<T> &)> initialization_function,
                std::function<double(std::vector<T> &)> evaluation_function,
                std::function<void(std::vector<T> &, double)> mutation_function,
                std::function<void(std::vector<T> &, std::vector<T> &)> crossover_function,
                unsigned int perfect_fitness) {

        // moved this here because to define initialization_function we need functionality from
        // population object, so we can't do it in constructor
        // TODO: think about a smarter solution
        initialize_genotypes(initialization_function);

        // Initial evaluation (this is called within the generation functions later)
        evaluate_all_fitnesses(evaluation_function);
        remember_best_genotype();

        for (unsigned int generation = 0; generation < max_generations; generation++) {
            create_new_population();
            crossover_population(crossover_function);
            mutate_population(mutation_function);
            report(generation);
            evaluate_all_fitnesses(evaluation_function);
            elitist();

            if (current_best_genotype.fitness == perfect_fitness) return;
        }
    }

    // Print the best variable values and the corresponding fitness
    void print_result() const {
        std::cout << "\n";
        std::cout << "  Best member:\n";
        std::cout << "\n";

        for (int i = 0; i < variable_count; i++) {
            std::cout << "  var(" << i << ") = " << current_best_genotype.genes[i] << "\n";
        }

        std::cout << "\n";
        std::cout << "  Best fitness = " << current_best_genotype.fitness << "\n";
    }

    int GetVariableCount() const {
        return variable_count;
    }

    Genotype<T> GetBestGenotype() const {
        return current_best_genotype;
    }
    // Return a random double between 0.0 and 1.0
    double GetRandomNormalizedDouble() {
        return random_normalized_double(mersenne_twister_engine);
    }

    // Return a random possible gene value at the given index, taking the lower and upper bounds
    // into account
    T GetRandomGeneValue(unsigned int index) {
        // using round because otherwise upperbound will only be returned when we get exactly 1.0
        // with GetRandomNormalizedDouble() when using int
        // TODO: does this work with doubles?
        return round(GetRandomNormalizedDouble() *
                         (upper_gene_bound[index] - lower_gene_bound[index]) +
                     lower_gene_bound[index]);
    }

  private:
    std::vector<Genotype<T>> genotypes;
    std::vector<Genotype<T>> new_genotypes;
    Genotype<T> current_best_genotype;

    /// Upper and lower bounds for the gene variables
    std::vector<T> upper_gene_bound;
    std::vector<T> lower_gene_bound;

    unsigned int variable_count;
    const unsigned int max_generations;
    const unsigned int size;

    // Probabilities
    const double crossover_probability;
    const double mutation_probability;

    // Random number generators
    std::random_device random_device;
    std::mt19937 mersenne_twister_engine;
    std::uniform_real_distribution<double> random_normalized_double;
    std::uniform_int_distribution<int> random_gene_index;

    // Parse a file with the format:
    // [variable_count]
    // [lower bound 1] [upper bound 1]
    // [lower bound 2] [upper bound 2]
    // ...
    // [lower bound variable_count] [upper bound variable_count]
    // for all variable_count variables.
    void parse_bounds_file(std::string filepath) {
        std::ifstream input;
        input.open(filepath);

        if (!input) {
            std::cerr << "Cannot open input file!" << std::endl;
            exit(1);
        }

        // reading in variable_count from file to make iniatilzation and testing easier
        input >> variable_count;
        lower_gene_bound.resize(variable_count);
        upper_gene_bound.resize(variable_count);
        // Parse and set bounds
        for (unsigned  int i = 0; i < variable_count; i++) {
            T lower_bound;
            T upper_bound;

            input >> lower_bound >> upper_bound;

            lower_gene_bound[i] = lower_bound;
            upper_gene_bound[i] = upper_bound;
        }

        input.close();
    }

    // Initialize Genotype genes to random values
    void initialize_genotypes(std::function<void(std::vector<T> &)> initialization_function) {
        for (Genotype<T> &genotype : genotypes) {
            genotype.set_variable_count(variable_count);
            initialization_function(genotype.genes);
        }
    }

    // Return a random int between 0 and the number of variables minus one, to be used as an index
    // for gene arrays
    unsigned int get_random_gene_index() {
        return random_gene_index(mersenne_twister_engine);
    }

    /// Iterates through the population and crosses over randomly selected pairs.
    void crossover_population(
        std::function<void(std::vector<T> &, std::vector<T> &)> crossover_function) {
        int current_match;
        int match_count = 0;

        for (unsigned int current_genotype = 0; current_genotype < size; ++current_genotype) {
            // Pick a random number for testing against the pre-defined crossover probability
            double x = GetRandomNormalizedDouble();

            if (x < crossover_probability) {
                // This genotype will be involved in a crossover
                ++match_count;

                if (match_count % 2 == 0) {
                    crossover(crossover_function, current_match, current_genotype);
                } else {
                    current_match = current_genotype;
                }
            }
        }
    }

    /// Performs a crossover on the two given parents.
    void crossover(std::function<void(std::vector<T> &, std::vector<T> &)> crossover_function,
                   int one, int two) {
        // Randomly select the point until which the crossover will be performed
        crossover_function(genotypes[one].genes, genotypes[two].genes);
    }

    // Elitist function: find the best and worse genotypes, remember the best, and replace the worst
    // by the previous best if the current best is worse than the previous best.
    void elitist() {
        int index_of_best;
        int index_of_worst;

        double best_fitness = genotypes[0].fitness;
        double worst_fitness = genotypes[0].fitness;

        // Find the best and the worst genotypes
        for (unsigned int i = 0; i < size - 1; ++i) {
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
    void evaluate_all_fitnesses(std::function<double(std::vector<T> &)> evaluation_function) {
        for (Genotype<T> &genotype : genotypes) {
            genotype.fitness = evaluation_function(genotype.genes);
        }
    }

    /// Returns the genotype with the highest fitness.
    void remember_best_genotype() {
        for (Genotype<T> &genotype : genotypes) {
            if (genotype.fitness > current_best_genotype.fitness) {
                // TODO: We might want to avoid the overhead of copying all the time?
                current_best_genotype = genotype;
            }
        }
    }

    // Mutate random genes of all genotypes.
    void mutate_population(std::function<void(std::vector<T> &, double)> mutation_function) {
        // Iterate over all genes in every genotype
        for (unsigned int i = 0; i < size; i++) {
            mutation_function(genotypes[i].genes, mutation_probability);
        }
    }

    // Selector function which generates a new population by keeping random good members of the
    // current one
    void create_new_population() {
        // Find the total fitness of the population
        double fitness_sum = 0.0;
        for (unsigned int i = 0; i < size; i++) {
            fitness_sum = fitness_sum + genotypes[i].fitness;
        }

        // Calculate the relative fitness of each member and save it
        for (unsigned int i = 0; i < size; i++) {
            genotypes[i].relative_fitness = genotypes[i].fitness / fitness_sum;
        }

        // Calculate the cumulative fitness of each member and save it
        genotypes[0].cumulative_fitness = genotypes[0].relative_fitness;
        for (unsigned int i = 1; i < size; i++) {
            genotypes[i].cumulative_fitness =
                genotypes[i - 1].cumulative_fitness + genotypes[i].relative_fitness;
        }

        // Select survivors using the cumulative fitnesses
        for (unsigned int i = 0; i < size; i++) {
            double p = GetRandomNormalizedDouble();

            if (genotypes[0].cumulative_fitness > p) {
                new_genotypes[i] = genotypes[0];
            } else {
                for (unsigned int j = 0; j < size - 1; j++) {
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
        for (unsigned int i = 0; i < size; i++) {
            genotypes[i] = new_genotypes[i];
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

        for (unsigned int i = 0; i < size; i++) {
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
};