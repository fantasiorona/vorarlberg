#include <vector>

/// A single member of a population
template <class T>
class Genotype {
  public:
    Genotype() = default;

    Genotype(unsigned int variable_count) {
        gene.resize(variable_count);
        upper.resize(variable_count);
        lower.resize(variable_count);
    }

    /// A set of variables
    std::vector<T> gene;

    /// Upper and lower bounds for the gene variables
    std::vector<T> upper;
    std::vector<T> lower;

    /// Fitness
    double fitness;
    double relative_fitness;
    double cumulative_fitness;
};