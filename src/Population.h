#include "Genotype.h"
#include <string>

template <class T>
class Population {
  public:
    std::vector<Genotype<T>> genotypes;
    std::vector<Genotype<T>> new_genotypes;

    unsigned int variable_count;
    unsigned int max_generations;
    unsigned int size;

    double crossover_probability;
    double mutation_probability;

    void crossover(int &seed);
    void elitist();
    void evaluate();
    int i4_uniform_ab(int a, int b, int &seed);
    void initialize(std::string filename, int &seed);
    void keep_the_best();
    void mutate(int &seed);
    double r8_uniform_ab(double a, double b, int &seed);
    void report(int generation);
    void selector(int &seed);
    void timestamp();
    void Xover(int one, int two, int &seed);
};