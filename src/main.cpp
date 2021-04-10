#include "Population.h"

#include <algorithm>
#include <chrono>
#include <map>

void print_magic_square(std::vector<int> &values);
void print_genes(std::vector<int> &values);
void create_input_data();
double evaluate_for_semi_magic_square(std::vector<int> &genes);
double evaluate_for_pandiagonal_magic_square(std::vector<int> &genes);
double evaluate_for_inlaid_semi_magic_square(std::vector<int> &genes);
void initialze_genotype(std::vector<int> &genes);
void mutate_genotype(std::vector<int> &genes, double mutationProbability);

// creating vectors here so we dont have to do it in each crossover function
std::vector<int> crossoverInverseSequence1;
std::vector<int> crossoverInverseSequence2;
std::vector<int> crossoverPositionSequence1;
std::vector<int> crossoverPositionSequence2;
void crossover_genotypes(std::vector<int> &, std::vector<int> &);
const int maxGenerations = 1000000000;
const int populationSize = 20;
const float crossoverChance = 0.8;
const float mutationChance = 0.05;
const std::string file = "inputs/normal_dimension_4.txt";
int dimension = 5;
// Create a new Population with genes consisting of three `double` variables
// global because we use it the mutation function
// yes this is spaghetti code
Population<int> population(file, maxGenerations, populationSize, crossoverChance, mutationChance);
int main(int argc, char **argv) {
    auto start = std::chrono::high_resolution_clock::now();

    // Uncomment to create input files
    // create_input_data();

    // Find dimension of magic square based on the input file
    int variableCount = population.GetVariableCount();
    dimension = sqrt(variableCount);

    // Check if input file makes sense for a square
    if (dimension * dimension != variableCount) {
        // TODO: Make this sentence make sense
        std::cerr << "ERROR: Input file doesnt contain a cleanly squared amount of values."
                  << std::endl;
        return -1;
    }
    // adjusting size according to variable count
    crossoverInverseSequence1.resize(variableCount);
    crossoverInverseSequence2.resize(variableCount);
    crossoverPositionSequence1.resize(variableCount);
    crossoverPositionSequence2.resize(variableCount);

    // Evolve the population to maximize the given equation
    if (std::stoi(argv[1]) == 1) {
        // Solve for normal semi-magic square
        population.evolve(initialze_genotype, evaluate_for_semi_magic_square, mutate_genotype,
                          crossover_genotypes, dimension * 2);
    } else if (std::stoi(argv[1]) == 2) {
        // Solve for pandiagonal magic square
        population.evolve(initialze_genotype, evaluate_for_pandiagonal_magic_square,
                          mutate_genotype, crossover_genotypes, dimension * 4);
    } else if (std::stoi(argv[1]) == 3) {
        // Solve for inlaid magic square
        population.evolve(initialze_genotype, evaluate_for_inlaid_semi_magic_square,
                          mutate_genotype, crossover_genotypes,
                          dimension * 2 + (dimension - 2) * 2);
    } else {
        std::cerr << "Unknown mode: " << argv[1] << std::endl;
    }
    // population.print_result();

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
    print_magic_square(bestMember.genes);
}
std::vector<int> calc_sums(const std::vector<int> &square, int dimension) {

    int total_sum = 0;
    std::vector<int> sums(dimension * 2 + 1);
    for (int row = 0; row < dimension; ++row) {
        for (int col = 0; col < dimension; ++col) {
            int val = square[row * dimension + col];
            sums[row] += val;
            sums[dimension + col] += val;
            total_sum += val;
        }
    }
    return sums;
}

std::vector<int> calc_cross_diagonals(const std::vector<int> &square, int dimension) {

    std::vector<int> sums(dimension * 2 + 1);
    for (int diag = 0; diag < dimension; ++diag) {
        for (int cell = 0; cell < dimension; ++cell) {
            int idx = diag + cell * (dimension - 1);
            if (diag + 1 <= cell) idx += dimension;
            sums[diag] += square[idx];

            idx = diag + cell * (dimension + 1);
            if (diag + cell >= dimension) idx -= dimension;
            sums[diag + dimension] += square[idx];
        }
    }
    return sums;
}

// this just uses the next smaller concentric square as it was the simplest to implement (and it
// already takes a lot of time to calculate)
std::vector<int> calc_inlaid_semi_magic_squares(const std::vector<int> &square, int dimension) {
    if (dimension <= 3) {
        return std::vector<int>(dimension * dimension, 0);
    }
    // create inline square from data
    int newDimension = (dimension - 2);

    std::vector<int> inlaid = std::vector<int>(newDimension * newDimension, 0);
    int idx = 0;
    for (int row = 1; row < dimension - 1; ++row) {
        for (int col = 1; col < dimension - 1; ++col) {
            inlaid[idx++] = square[row * dimension + col];
        }
    }

    // passing a new vector is pretty sub optimal
    // TODO: pass original vector with a offset value
    return calc_sums(inlaid, newDimension); // semi-magic squares
}

int count_correct_sums(const std::vector<int> &sums, int dimension) {
    // from:
    // https://www.dr-mikes-math-games-for-kids.com/magic-square-magic-number.html#:~:text=A%20magic%20square%20is%20a,number"%20of%20the%20magic%20square
    int magicConstant = (dimension * dimension + 1) * dimension / 2;
    int correctSums = 0;
    for (auto sum : sums) {
        if (sum == magicConstant) {
            ++correctSums;
        }
    }
    return correctSums;
}

double evaluate_for_semi_magic_square(std::vector<int> &genes) {
    std::vector<int> semiMagicSums = calc_sums(genes, dimension);
    int correctSums = count_correct_sums(semiMagicSums, dimension);
    return correctSums;
}

double evaluate_for_pandiagonal_magic_square(std::vector<int> &genes) {

    std::vector<int> semiMagicSums = calc_sums(genes, dimension);

    std::vector<int> panDiagonalSums = calc_cross_diagonals(genes, dimension);

    int correctSums = 0;
    correctSums += count_correct_sums(semiMagicSums, dimension);
    correctSums += count_correct_sums(panDiagonalSums, dimension);

    return correctSums;
}

double evaluate_for_inlaid_semi_magic_square(std::vector<int> &genes) {

    std::vector<int> semiMagicSums = calc_sums(genes, dimension);
    std::vector<int> inlaidSums = calc_inlaid_semi_magic_squares(genes, dimension);

    int correctSums = 0;
    correctSums += count_correct_sums(semiMagicSums, dimension);

    // we cant simply compare inlaid square with magic constant because the values are not only
    // normal (meaning for example 1-9 for dimension 3), but can contain values from the bigger
    // encircling square (so 1-25 for dimension 3)
    // instead we count the duplicate values in the sums to find how good the inlaid square is
    // from https://stackoverflow.com/a/32238115
    std::map<int, int> sumDuplicates;
    for_each(inlaidSums.begin(), inlaidSums.end(),
             [&sumDuplicates](int val) { sumDuplicates[val]++; });
    // finding max duplicates
    int sumDuplicateMax = 0;
    for (auto duplicate : sumDuplicates) {
        int count = duplicate.second;
        sumDuplicateMax = sumDuplicateMax < count ? count : sumDuplicateMax;
    }
    correctSums += sumDuplicateMax;

    return correctSums;
}

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

// This only handles normal magic square
// TODO: also implement for non normal magic squares
void mutate_genotype(std::vector<int> &genes, double mutationProbability) {
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

// Adapted from https://github.com/guisehn/genetic-magic-square-finder
// TODO: maybe put inverse and position vectors in a twodimensional vector so we reduce code
// duplication
void crossover_genotypes(std::vector<int> &genes1, std::vector<int> &genes2) {
    // calculate inversion sequences
    int geneAmount = genes1.size();
    // inverse sequence counts the amount of bigger values to the left of the current value and
    // stores it at the corresponding index
    // so for example for the gene [2 3 4 5 1 6 7 8 9] we start with the value 1
    // it is at index 5 and there are 4 bigger values to the left
    // so the value in the inverse sequence at index 1 would be 4
    // this uniquely identifies the gene but also allows to cross with out producing duplicates
    // because translating back a inverse sequence only allows the creation of different values for
    // each gene variable

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

void print_genes(std::vector<int> &values) {
    for (unsigned int i = 0; i < values.size(); ++i) {
        if (values[i] < 0) {
            std::cout << "X, ";
        } else {
            std::cout << values[i] << ", ";
        }
    }
    std::cout << std::endl;
}

// do not look at this function it is one of the worst things i have ever written
void print_magic_square(std::vector<int> &values) {
    std::cout << std::endl << "Magic Square:" << std::endl << std::endl;
    int magicConstant = (dimension * dimension + 1) * dimension / 2;
    std::cout << "Magic constant should be: " << magicConstant << std::endl;

    std::vector<int> inlaidSums = calc_inlaid_semi_magic_squares(values, dimension);
    std::map<int, int> sumDuplicates;
    for_each(inlaidSums.begin(), inlaidSums.end(),
             [&sumDuplicates](int val) { sumDuplicates[val]++; });
    // finding max duplicates
    int sumDuplicateMax = 0;
    int bestMagicConstant = -1;
    for (auto duplicate : sumDuplicates) {
        int count = duplicate.second;
        if (sumDuplicateMax < count) {
            sumDuplicateMax = count;
            bestMagicConstant = duplicate.first;
        }
    }
    int inlaidCorrectSums = 0;
    for (auto s : inlaidSums) {
        if (s == bestMagicConstant) ++inlaidCorrectSums;
    }
    std::cout << "Magic constant with most duplicates in the inlaid square is: "
              << bestMagicConstant << std::endl
              << std::endl;

    int semiMagicCorrectSums = 0;
    for (unsigned int row = 0; row < values.size(); row += dimension) {
        int sum = 0;
        for (int col = 0; col < dimension; ++col) {
            int curValue = values[row + col];
            if (curValue < 100) {
                std::cout << "0";
            }
            if (curValue < 10) {
                std::cout << "0";
            }
            std::cout << curValue << " ";
            sum += curValue;
        }

        std::cout << "- ";
        if (sum < 100) {
            std::cout << "0";
        }
        if (sum < 10) {
            std::cout << "0";
        }
        std::cout << sum << std::endl;
        if (sum == magicConstant) {
            ++semiMagicCorrectSums;
        }
    }
    std::cout << " ";
    for (int i = 0; i < dimension; i++) {
        std::cout << "|   ";
    }
    std::cout << std::endl;
    for (int col = 0; col < dimension; ++col) {
        int sum = 0;
        for (unsigned int row = 0; row < values.size(); row += dimension) {
            int curValue = values[row + col];
            sum += curValue;
        }
        if (sum < 100) {
            std::cout << "0";
        }
        if (sum < 10) {
            std::cout << "0";
        }
        std::cout << sum << " ";
        if (sum == magicConstant) {
            ++semiMagicCorrectSums;
        }
    }

    std::cout << std::endl << std::endl;
    std::cout << "Broken diagonal sums are: " << std::endl;
    std::vector<int> pandiagonalSums = calc_cross_diagonals(values, dimension);
    std::cout << "Left to Right: ";
    for (int i = dimension; i < dimension * 2; ++i) {
        std::cout << std::to_string(pandiagonalSums[i]) << " ";
    }
    std::cout << std::endl;
    std::cout << "Right to Left: ";
    for (int i = 0; i < dimension; ++i) {
        std::cout << std::to_string(pandiagonalSums[i]) << " ";
    }
    int pandiagonalCorrectSums = count_correct_sums(pandiagonalSums, dimension);

    std::cout << std::endl << std::endl;
    std::cout << "Inlaid sums are: " << std::endl;
    std::cout << "Rows: ";
    for (int i = 0; i < dimension - 2; i++) {
        std::cout << inlaidSums[i] << ", ";
    }
    std::cout << std::endl;
    std::cout << "Columns: ";
    for (int i = dimension - 2; i < (dimension - 2) * 2; i++) {
        std::cout << inlaidSums[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << std::endl << std::endl;
    std::cout << "Amount of correct semi magic sums: " << semiMagicCorrectSums << std::endl;
    std::cout << "Amount for semi magic square should be: " << std::to_string(dimension * 2)
              << std::endl
              << std::endl;
    std::cout << "Amount of correct pandiagonal sums: " << semiMagicCorrectSums << "(semi) + "
              << pandiagonalCorrectSums
              << "(pandiagonal) = " << semiMagicCorrectSums + pandiagonalCorrectSums << std::endl;
    std::cout << "Amount for pandiagonal magic square should be: " << std::to_string(dimension * 4)
              << std::endl
              << std::endl;
    std::cout << "Amount of correct inlaids sums: " << semiMagicCorrectSums << "(semi) + "
              << inlaidCorrectSums << "(inlaid) = " << semiMagicCorrectSums + inlaidCorrectSums
              << std::endl;
    std::cout << "Amount for inlaids magic square should be: "
              << std::to_string(dimension * 2 + (dimension - 2) * 2) << std::endl
              << std::endl;
}

void create_input_data() {
    // File writing from
    // https://www.sqlnethub.com/blog/how-to-write-to-text-file-from-c-plus-plus-program/

    // create data for each necessary dimension
    for (int dim = 3; dim <= 9; ++dim) {
        // maxvalue for a normal dataset is dimension sqaured
        int maxValue = dim * dim;
        // creating string for lower and upper bounds in necessary format
        std::string data = "1 " + std::to_string(maxValue);
        // open file for writing
        std::ofstream fw(("inputs\\normal_dimension_" + std::to_string(dim) + ".txt"),
                         std::ofstream::out);
        if (fw.is_open()) {
            fw << std::to_string(maxValue) << std::endl;
            // writing bounds as a line for every gene that is later neede
            for (int line = 0; line < maxValue; ++line) {
                fw << data << std::endl;
            }
            fw.close();
        } else
            std::cout << "Problem with opening file";
        for (int factor = 2; factor < 5; ++factor) {
            int nonNormalMaxvalue = maxValue * factor;
            std::string nonNormalData = "1 " + std::to_string(nonNormalMaxvalue);
            // open file for writing
            std::ofstream fw(("inputs\\notnormal_factor_" + std::to_string(factor) + "_dimension_" +
                              std::to_string(dim) + ".txt"),
                             std::ofstream::out);
            if (fw.is_open()) {
                fw << std::to_string(maxValue) << std::endl;
                for (int j = 0; j < maxValue; j++) {
                    fw << nonNormalData << "\n";
                }
                fw.close();
            } else
                std::cout << "Problem with opening file";
        }
    }
}