#include <limits>
#include <vector>

/// A single member of a population
template <class T>
class Genotype {
  public:
    Genotype()
        : fitness(-std::numeric_limits<float>::infinity()), relative_fitness(0.0),
          cumulative_fitness(0.0) {
    }

    void set_variable_count(unsigned int variable_count) {
        genes.resize(variable_count);
    }

    Genotype<T>& operator=(const Genotype<T>& other) {
      fitness = other.fitness;
      relative_fitness = other.relative_fitness;
      cumulative_fitness = other.cumulative_fitness;
      genes = other.genes;

      return *this;
    }

    /// A set of variables
    std::vector<T> genes;

    /// Fitness
    double fitness;
    double relative_fitness;
    double cumulative_fitness;
};
