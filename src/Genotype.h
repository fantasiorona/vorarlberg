#include <vector>

/// A single member of a population
template <class T>
class Genotype {
  public:
    Genotype() : fitness(0.0), relative_fitness(0.0), cumulative_fitness(0.0) {
    }

    void set_variable_count(unsigned int variable_count) {
        genes.resize(variable_count);
    }

    /// A set of variables
    std::vector<T> genes;

    /// Fitness
    double fitness;
    double relative_fitness;
    double cumulative_fitness;
};