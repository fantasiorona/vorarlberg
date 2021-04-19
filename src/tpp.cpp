/* ==============================================================================================
TPP
Implement a GP version of the Traveling Salesman Problem given as a Traveling Purchaser Problem.

where:
- each city has a market with just one good
- this good has a price equal to the length of the name of the city.
- Solve the problem for the Romanian roads given.

# references
- https://github.com/weilheim/Genetic-Algorithm-for-TSP
- https://www.geeksforgeeks.org/traveling-salesman-problem-using-genetic-algorithm/

# papers
- DESIGN AND IMPLEMENTATION OF A PARALLEL GENETIC ALGORITHM FOR THE TRAVELLING PURCHASER PROBLEM -
Luiz S. Ochi, Lficia M. A. Drummond
- The Generalized Traveling Salesman Problem: A New Genetic Algorithm Approach - John Silberholz and
Bruce L. Golden
- A NEW GENETIC ALGORITHM APPLIED TO THE TRAVELING SALESMAN PROBLEM - Sawsan K. Amous, Taicir
Loukil, Semya Elaoud, Clarisse Dhaenens
- The Traveling Purchaser Problem and its variants - Daniele Manerbaa, Renata Mansini, Jorge
Riera-Ledesma
- Transgenetic algorithm for the Traveling Purchaser Problem- M.C. Goldbarg, L.B. Bagi, E.F.G.
Goldbarg

# problems:
- crossover: how to handle different array sizes?
    --> building clusters which connect all purchasable goods, so the size is always = provided
goods
- mutation: how to quarantee that swapped city is reachable?
    --> clustering also solves this because distances between all traversable cities are stored
- need to create a "circle" with a given starting point, returning to after collecting all goods?
    --> implemented using -c cmd argument
- necessary to provide a starting point or just start from city with first good to purchase?
    --> allow providing a start city cmd argument with --from <idx> or <cityname>



# examples:

cmd: ./tpp -i 100 -p 10 1 2 3 4 5 6 7 8 9 10 11 12 13 13

provided:           14, 12, 19, 18, 11, 10, 4,

12 - 14 = 3086: bucharest - urziceni: 85
12 - 19 = 3091: bucharest - neamt: 406
18 - 19 = 4627: iasi - neamt: 87
11 - 18 = 2834: pitesti - iasi: 420
10 - 11 = 2571: rimnicuvilcea - pitesti: 97
4 - 10 = 1034: oradea - rimnicuvilcea: 231
------------------------------------------
1448

another case:       0, 1, 2, 10, 12, 14, 15,

0 - 1 = 1: arad - zerind: 75
1 - 2 = 258: zerind - sibiu: 215
2 - 10 = 522: sibiu - rimnicuvilcea: 80
10 - 12 = 2572: rimnicuvilcea - bucharest: 198
12 - 14 = 3086: bucharest - urziceni: 85
14 - 15 = 3599: urziceni - hirsova: 98
------------------------------------------
751


min found with 1 iteration an 10000 population:
for genome:         2, 10, 11, 12, 14, 17, 18

2 - 10 = 522: sibiu - rimnicuvilcea: 80
10 - 11 = 2571: rimnicuvilcea - pitesti: 97
11 - 12 = 2828: pitesti - bucharest: 101
12 - 14 = 3086: bucharest - urziceni: 85
14 - 17 = 3601: urziceni - vaslui: 142
17 - 18 = 4370: vaslui - iasi: 92
------------------------------------------
597

==============================================================================================*/

#include "Population.h"
#include "utility.h"

#include <cstring>
#include <functional>
#include <set>
#include <unordered_set>

#define VERBOSE
//#define FILE_OUT

// using key for distance map as bitmask of idxFrom and idxTo
// e.g. from node 6 to 10 is in binary 0110 and 1010 which results in an unique key of
// 0000 0110    0000 1010 -> 1546 in dec
#define KEY_UPPER 8    // 8 bitshift is enough for 255 nodes, we only have 20
#define KEY_LOWER 0xFF // pow(2, KEY_UPPER)

typedef size_t VisitedCity;
typedef std::vector<VisitedCity> CitySequence;

void print_sequence(CitySequence &values) {
    // TODO: Consider printing names rather than indices?
    for (size_t value : values) {
        std::cout << value << ", ";
    }
    std::cout << std::endl;
}

std::function<double(CitySequence &)> get_evaluation_function(std::map<int, int> &city_distances,
                                                              bool circle) {
    return [city_distances, circle](CitySequence &genes) {
        int sum = 0, from, to;

        std::vector<size_t>::iterator it_to;
        for (auto it_from = genes.begin(); it_from != genes.end(); ++it_from) {
            from = *it_from;
            it_to = std::next(it_from, 1);
            if (it_to == genes.end()) {
                if (false == circle) break;

                // if purchaser has to travel back to start add closest distance
                it_to = genes.begin();
            }

            to = *it_to;
            if (from > to)
                std::swap(from, to); // assure ascending direction to find in distance map

            sum += city_distances.at(from << KEY_UPPER | to);
        }

        // Return the negative sum because a higher fitness is considered better --> we're
        // optimizing towards 0
        return -sum;
    };
}

std::function<void(CitySequence &)>
get_initialization_function(std::map<int, std::vector<size_t>> cluster, int idx_first_city) {
    return [cluster, idx_first_city](CitySequence &genes) {
        for (const auto &entry : cluster) {
            // choose random index for multiple options
            int idx = (entry.second.size() == 1) ? 0 : rand() % entry.second.size();
            // genome.insert(entry.second[idx]);
            genes.push_back(entry.second[idx]);
        }

        std::random_shuffle(genes.begin() + (-1 != idx_first_city),
                            genes.end()); // don't shuffle starting point if given
    };
}

// Declared here so they're re-used in each crossover sequenc
std::vector<int> crossover_inverse_sequence1;
std::vector<int> crossover_inverse_sequence2;
std::vector<int> crossover_position_sequence1;
std::vector<int> crossover_position_sequence2;
void crossover_genotypes(CitySequence &genes1, CitySequence &genes2,
                         Population<VisitedCity> &population) {
    // calculate inversion sequences
    size_t gene_amount = genes1.size();
    // inverse sequence counts the amount of bigger values to the left of the current value and
    // stores it at the corresponding index
    // so for example for the gene [2 3 4 5 1 6 7 8 9] we start with the value 1
    // it is at index 5 and there are 4 bigger values to the left
    // so the value in the inverse sequence at index 1 would be 4
    // this uniquely identifies the gene but also allows to cross with out producing duplicates
    // because translating back a inverse sequence only allows the creation of different values for
    // each gene variable

    // iterate over each value
    for (size_t i = 0; i < gene_amount; ++i) {
        // initialize to zero
        crossover_inverse_sequence1[i] = 0;
        crossover_inverse_sequence2[i] = 0;
        bool found1 = false;
        bool found2 = false;
        // for each value we start at the left
        for (size_t j = 0; j < gene_amount; ++j) {
            // and continue until we find the value
            if (!found1) {
                // if we pass a bigger value we increment
                if (genes1[j] > i + 1) {
                    crossover_inverse_sequence1[i] = crossover_inverse_sequence1[i] + 1;
                }
                // stop when we find the value
                if (genes1[j] == i + 1) {
                    found1 = true;
                }
            }
            // do same for the second parent
            if (!found2) {
                if (genes2[j] > i + 1) {
                    crossover_inverse_sequence2[i] = crossover_inverse_sequence2[i] + 1;
                }
                if (genes2[j] == i + 1) {
                    found2 = true;
                }
            }
        }
    }

    // crossover afterwards is a simple swap operation with a randomly defined cutoff point
    size_t cutoff_point = population.GetRandomGeneValue(0);
    for (size_t i = cutoff_point; i < gene_amount; ++i) {
        std::swap(crossover_inverse_sequence1[i], crossover_inverse_sequence2[i]);
    }
    // translate back to actual genes
    // we start with creating a position vector
    // this happens by counting greater or equal elements to the right of the current value
    // this is basically the inverse operation to the inverse sequence creation and gives the
    // information where each value is positioned in the gene vector
    for (int i = gene_amount - 1; i >= 0; --i) {
        // initialaze to current inverse value
        crossover_position_sequence1[i] = crossover_inverse_sequence1[i];
        crossover_position_sequence2[i] = crossover_inverse_sequence2[i];
        // for each value to the right
        for (size_t j = i + 1; j < gene_amount; ++j) {
            // increase if its bigger or equal
            if (crossover_position_sequence1[j] >= crossover_position_sequence1[i]) {
                crossover_position_sequence1[j] = crossover_position_sequence1[j] + 1;
            }
            // do the same for the other parent
            if (crossover_position_sequence2[j] >= crossover_position_sequence2[i]) {
                crossover_position_sequence2[j] = crossover_position_sequence2[j] + 1;
            }
        }
    }

    // actually create genes by using the position information to find the position of the value
    // add 1 because values are from 1-9 while vector indizes are from 0-8
    for (size_t i = 0; i < gene_amount; ++i) {
        genes1[crossover_position_sequence1[i]] = i + 1;
        genes2[crossover_position_sequence2[i]] = i + 1;
    }
}

void mutate_genotype(CitySequence &genes, double mutation_probability,
                     Population<VisitedCity> &population) {
    int variable_amount = genes.size();

    // iterating over each gene
    for (int j = 0; j < variable_amount; ++j) {
        // random chance to mutate gene
        double x = population.GetRandomNormalizedDouble();

        if (x < mutation_probability) {
            int swapIndex = population.GetRandomGeneValue(0);
            // make sure we dont swap in place
            while (swapIndex == j) {
                ++swapIndex;
                if (swapIndex == variable_amount + 1) swapIndex = 1;
            }
            std::swap(genes[j], genes[swapIndex]);
        }
    }
}

int main(int argc, char *argv[]) {
    size_t iterations = 1;      // 200;    // extends the population for every iteration?
    size_t populationSize = 20; // how many different routes are used
    int mutationRate = 5;       // percentage of which a mutation should happen
    int crossoverRate = 70;     // percentage of how often crossover should happen
    bool circle = false;        // flag wheter traveling from last good back to first should happen
    std::set<int> goods;        // provided goods to buy on tour
    int startIdx = -1; // starting point index for route, if not provided just start from first good
    std::string startName = ""; // starting point cityname
    std::string dataPath = "inputs/romaniaroads.pl";

    // 1. read input params
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    for (int i = 1; i < argc; ++i) // first arg is app name
    {
        if (strcmp(argv[i], "-f") == 0)
            dataPath = argv[++i];
        else if (strcmp(argv[i], "-i") == 0)
            iterations = std::stoi(argv[++i]);
        else if (strcmp(argv[i], "-p") == 0)
            populationSize = std::stoi(argv[++i]);
        else if (strcmp(argv[i], "-c") == 0)
            circle = true;
        else if (strcmp(argv[i], "-m") == 0)
            mutationRate = std::stoi(argv[++i]);
        else if (strcmp(argv[i], "-x") == 0)
            crossoverRate = std::stoi(argv[++i]);
        else if (strcmp(argv[i], "--from") == 0) {
            try {
                startIdx = std::stoi(argv[++i]);
            } catch (const std::invalid_argument &) {
                startName = argv[i];
            }
        } else
            goods.insert(std::stoi(argv[i])); // remaining given values are interpreted as goods
                                              // goods are names for cities -> using indices
    }

    // 2. read given data containing cities and streets
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    // const std::vector<city> &cities = readInput("data/romaniaroads.pl");
    readInput(dataPath);

    size_t idx = 0;
#ifdef VERBOSE
    std::cout << "cities with connected roads (" << cities.size()
              << "), index is considered as provided good: " << std::endl;
    for (const auto &city : cities) {
        std::cout << "  " << idx++ << ": " << city.name << " (cost: " << city.good
                  << ") with roads (" << city.roads.size() << "): " << std::endl;
        for (Road road : city.roads) {
            std::cout << "\t" << cities[road.to].name << ": " << road.distance << std::endl;
        }
    }
#endif

    if (cities.size() < 1) {
        std::cout << "no valid city network found, please contact your local administrator."
                  << std::endl;
        return 0;
    }

    while (goods.size() < 2) {
        std::cout << "not enough goods provided to justify some work: " << goods.size()
                  << " min: 2, given: " << goods.size() << std::endl
                  << "please provide indices for cities (0 - " << (cities.size() - 1)
                  << "): <idx1> <idx2> ..." << std::endl;

        std::string str = "";
        std::getline(std::cin, str);
        std::vector<std::string> input = split(str, " ");
        for (const std::string &idx_city : input) {
            goods.insert(std::stoi(idx_city));
        }
    }

    // select cluster nodes which provide given goods
    idx = 0;
    std::map<int, std::vector<size_t>> cluster;
    // for (const auto &city : cities) {
    for (size_t idx_city = 0; idx_city < cities.size(); ++idx_city) {
        // if (std::find(goods.begin(), goods.end(), city.good) != goods.end()) {
        if (std::find(goods.begin(), goods.end(), idx_city) != goods.end()) {
            if (cluster.count(idx_city)) {
                cluster[idx_city].push_back(idx);
            } else {
                cluster.insert(std::make_pair(idx_city, std::vector<size_t>{idx}));
            }
        }
        ++idx;
    }

    // check for starting node -> prepend to cluster if provided
    if (false == startName.empty()) {
        if (cityOrdinals.count(startName)) {
            startIdx = cityOrdinals[startName];
        } else {
            std::cout << "given starting point \"" << startName << "\" invalid!" << std::endl;
            return 0;
        }
    }

    if (-1 != startIdx) {
        if (cities.size() > (size_t)startIdx) {
            cluster.insert(std::make_pair(0, std::vector<size_t>{(size_t)startIdx}));
            // remove starting point good from cluster because is already collected
            if (cluster.count(cities[startIdx].good)) {
                cluster.erase(cities[startIdx].good);
            }
        } else {
            std::cout << "starting point index \"" << startIdx << "\" invalid!" << std::endl;
            return 0;
        }
    }

    if (cluster.size() < 2) {
        std::cout << "not enough valid goods provided to justify some work: " << cluster.size()
                  << " min: 2" << std::endl;
        return 0;
    }

#ifdef VERBOSE
    std::cout << "selected cluster nodes (" << cluster.size() << "):" << std::endl;
    for (const auto &entry : cluster) {
        std::cout << "  good " << entry.first << ":" << std::endl;
        for (int node : entry.second) {
            std::cout << "\t" << node << " - " << cities[node].name << " (" << cities[node].good
                      << ")" << std::endl;
        }
    }
#endif

    // calculate distance for cluster nodes (could be done hardcoded because we only use the
    // romaniaroads data file)
    std::map<int, int> node_distances; // key is bitwise combination of node indices, only storing
                                       // one direction (smaller index to bigger)
    // omg noooe what am i doin? oO
    // TODO: double iterate cities, starting inner loop at current outer loop index, using only
    // relevant nodes and discarding same goods
    for (const auto &setFrom : cluster) {
        for (const auto &setTo : cluster) {
            if (setFrom.first == setTo.first)
                continue; // no connections needed between cities with the same goods

            for (int from : setFrom.second) {
                for (int to : setTo.second) {
                    if (from >= to) continue; // ignore duplicates and reverse order

                    int key = (from << KEY_UPPER | to);
                    int distance = minRoute(&cities[from], &cities[to]);
                    node_distances.insert(std::make_pair(key, distance));
                }
            }
        }
    }

#ifdef VERBOSE
    std::cout << "calculated distances (" << node_distances.size()
              << "):" << std::endl; // should be "factorial for sums" - cities with the same good
    for (const auto &dist : node_distances) {
        std::cout << "  " << (dist.first >> KEY_UPPER) << " - " << (dist.first & KEY_LOWER) << " = "
                  << dist.first << ": " << cities[dist.first >> KEY_UPPER].name << " - "
                  << cities[dist.first & KEY_LOWER].name << ": " << dist.second << std::endl;
    }
#endif

    // Get list of allowed values (actually not really necessary because we only swap?)
    std::vector<VisitedCity> allowed_values;
    for (const auto &entry : cluster) {
        allowed_values.emplace_back(entry.first);
    }

    // adjusting size according to variable count
    crossover_inverse_sequence1.resize(allowed_values.size());
    crossover_inverse_sequence2.resize(allowed_values.size());
    crossover_position_sequence1.resize(allowed_values.size());
    crossover_position_sequence2.resize(allowed_values.size());

    // FIXME: It'd be unnecessary to pass a file with min/max values here. Give the option of
    // initializing min/max to set values for all Genotypes!
    Population<VisitedCity> population(allowed_values, 500, populationSize,
                                       static_cast<float>(crossoverRate) / 100.0,
                                       static_cast<float>(mutationRate) / 100.0);

    // Fitnesses are negative, so 0 (meaning no distance is traversed at all) would be ideal
    population.evolve(get_initialization_function(cluster, startIdx),
                      get_evaluation_function(node_distances, circle), mutate_genotype,
                      crossover_genotypes, 0);

    return 0;
}
