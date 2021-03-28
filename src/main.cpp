#include "Genotype.h"
#include "Population.h"
#include <fstream>
#include <map>
#include <vector>

void print_magic_square(std::vector<int> values);
void create_input_data();
double evaluate_for_semi_magic_square(std::vector<int> genes);
double evaluate_for_pandiagonal_magic_square(std::vector<int> genes);

void mutate_genotype(std::vector<int> &genes, double mutationProbability);
int dimension = 0;
// Create a new Population with genes consisting of three `double` variables
// global because we use it the mutation function
// yes this is spaghetti code
Population<int> population("inputs\\normal_dimension_4.txt", 100, 2500, 0.4, 0.8);
int main() {

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

    // Evolve the population to maximize the given equation
    population.evolve(evaluate_for_semi_magic_square, mutate_genotype);

    // Print the result
    population.print_result();
    Genotype<int> bestMember = population.GetBestGenotype();
    print_magic_square(bestMember.genes);
}

std::vector<int> calc_sums(const std::vector<int> &square, int dimension) {

    int total_sum = 0;
    std::vector<int> sums(dimension * 2);
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

    std::vector<int> sums(dimension * 2); // /

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

double evaluate_for_semi_magic_square(std::vector<int> genes) {
    std::vector<int> semiMagicSums = calc_sums(genes, dimension);
    int correctSums = count_correct_sums(semiMagicSums, dimension);
    return correctSums;
}

double evaluate_for_pandiagonal_magic_square(std::vector<int> genes) {

    std::vector<int> semiMagicSums = calc_sums(genes, dimension);

    std::vector<int> panDiagonalSums = calc_cross_diagonals(genes, dimension);

    int correctSums = 0;
    correctSums += count_correct_sums(semiMagicSums, dimension);
    correctSums += count_correct_sums(panDiagonalSums, dimension);

    return correctSums;
}

// This only handels normal magic square
// TODO: also implement for non normal magic squares
void mutate_genotype(std::vector<int> &genes, double mutationProbability) {
    int variableCount = genes.size();

    // TODO: creating this vectors here probably isn't ideal
    std::vector<bool> usedValues;
    for (int i = 0; i < variableCount; ++i) {
        usedValues.push_back(false);
    }

    // iterating over each gene
    for (int j = 0; j < variableCount; ++j) {
        // random chance to mutate gene
        double x = population.GetRandomNormalizedDouble();
        int newValue = 0;
        if (x < mutationProbability) {
            // get a random value
            newValue = population.GetRandomGeneValue(j);
        } else {
            // take the old value when not mutating
            newValue = genes[j];
        }

        // count up while its not valid
        while (usedValues[newValue - 1]) {
            ++newValue;
            newValue = (newValue == variableCount + 1 ? 1 : newValue);
        }
        // switching used flag
        usedValues[newValue - 1] = true;
        genes[j] = newValue;
    }
}

void print_magic_square(std::vector<int> values) {
    std::cout << std::endl << "Magic Square:" << std::endl << std::endl;
    int magicConstant = (dimension * dimension + 1) * dimension / 2;
    std::cout << "Magic constant should be: " << magicConstant << std::endl << std::endl;
    int semiMagicCorrectSums = 0;
    int pandiagonalCorrectSums = 0;
    for (int row = 0; row < values.size(); row += dimension) {
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
        for (int row = 0; row < values.size(); row += dimension) {
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

    // TODO: code duplication from evaluate_for_pandiagonal_magic_square()
    std::vector<int> leftToRightSums;
    for (int col = 0; col < dimension; ++col) {
        int sum = 0;
        int offset = 0;
        for (int step = 0; step < dimension; ++step) {
            sum += values[col + offset];
            if (dimension - 1 - col == step) {
                offset += 1;
            } else {
                offset += (dimension + 1);
            }
        }
        leftToRightSums.push_back(sum);
        if (sum == magicConstant) {
            ++pandiagonalCorrectSums;
        }
    }
    std::vector<int> rightToLeftSums;
    for (int col = dimension - 1; col >= 0; --col) {
        int sum = 0;
        int offset = 0;
        for (int step = 0; step < dimension; ++step) {
            sum += values[col + offset];
            if (col == step) {
                offset += (2 * dimension - 1);
            } else {
                offset += (dimension - 1);
            }
        }
        rightToLeftSums.push_back(sum);
        if (sum == magicConstant) {
            ++pandiagonalCorrectSums;
        }
    }

    std::cout << "Left to Right: ";
    for (auto s : leftToRightSums) {
        std::cout << std::to_string(s) << " ";
    }
    std::cout << std::endl;
    std::cout << "Right to Left: ";
    for (auto s : rightToLeftSums) {
        std::cout << std::to_string(s) << " ";
    }
    std::cout << std::endl << std::endl;
    std::cout << "Amount of correct semi magic sums: " << semiMagicCorrectSums << std::endl;
    std::cout << "Amount for semi magic square should be: " << std::to_string(dimension * 2)
              << std::endl;
    std::cout << "Amount of correct pandiagonal sums: "
              << semiMagicCorrectSums + pandiagonalCorrectSums << std::endl;
    std::cout << "Amount for pandiagonal magic square should be: " << std::to_string(dimension * 4);
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