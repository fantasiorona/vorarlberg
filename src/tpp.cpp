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

            // FIXME: Debug print until the crossover function is fixed (visualizes invalid values)
            try {
                sum += city_distances.at(from << KEY_UPPER | to);
            } catch (std::exception &e) {
                std::cout << "error: invalid distance from " << from << " to " << to << std::endl;
                throw e;
            }
        }

        // Return the negative sum because a higher fitness is considered better --> we're
        // optimizing towards 0
        return -sum;
    };
}

std::function<void(CitySequence &)>
get_initialization_function(std::set<VisitedCity> &allowed_values, int idx_first_city) {
    return [allowed_values, idx_first_city](CitySequence &genes) {
        // TODO: better use std::copy?
        size_t idx = 0;
        for (const VisitedCity &node : allowed_values) {
            genes[idx] = node;
            if (node == (size_t)idx_first_city) std::swap(genes[0], genes[idx]);
            ++idx;
        }
        std::random_shuffle(genes.begin() + (-1 != idx_first_city),
                            genes.end()); // don't shuffle starting point if given
#ifdef VERBOSE
        // std::cout << "initial population: " << std::endl;
        // print_sequence(genes);
#endif
    };
}

// Declared here so they're re-used in each crossover sequence
std::vector<int> crossover_inverse_sequence1;
std::vector<int> crossover_inverse_sequence2;
std::vector<int> crossover_position_sequence1;
std::vector<int> crossover_position_sequence2;
std::vector<size_t> gene_indices1;
std::vector<size_t> gene_indices2;
std::function<void(CitySequence &, CitySequence &, Population<VisitedCity> &)>
get_crossover_function(std::set<VisitedCity> &allowed_values) {

    return [allowed_values](CitySequence &genes1, CitySequence &genes2,
                            Population<VisitedCity> &population) {
        size_t gene_amount = genes1.size();
        // transform genes into a form where we have continous values (meaning from 0 to
        // gene_amount-1), because otherwise this form of crossover does not make that much sense
        // as it transform the original values into a form where each value is defined by its
        // position in the vector. This only works if we know exactly which values have been in the
        // vector before hand, otherwise the inversion sequences are not unambiguous.
        for (size_t i = 0; i < gene_amount; i++) {
            gene_indices1[i] =
                std::distance(allowed_values.begin(), allowed_values.find(genes1[i]));
            gene_indices2[i] =
                std::distance(allowed_values.begin(), allowed_values.find(genes2[i]));
        }
        // calculate inversion sequences
        // inverse sequence counts the amount of bigger values to the left of the current value and
        // stores it at the corresponding index
        // so for example for the gene [2 3 4 5 1 6 7 8 9] we start with the value 1
        // it is at index 5 and there are 4 bigger values to the left
        // so the value in the inverse sequence at index 1 would be 4
        // this uniquely identifies the gene but also allows to cross with out producing duplicates
        // because translating back a inverse sequence only allows the creation of different values
        // for each gene variable

        // iterate over each value
        for (int i = 0; i < gene_amount; ++i) {
            // initialize to zero
            crossover_inverse_sequence1[i] = 0;
            crossover_inverse_sequence2[i] = 0;
            bool found1 = false;
            bool found2 = false;
            // for each value we start at the left
            for (int j = 0; j < gene_amount; ++j) {
                // and continue until we find the value
                if (!found1) {
                    // if we pass a bigger value we increment
                    if (gene_indices1[j] > i) {
                        crossover_inverse_sequence1[i] = crossover_inverse_sequence1[i] + 1;
                    }
                    // stop when we find the value
                    if (gene_indices1[j] == i) {
                        found1 = true;
                    }
                }
                // do same for the second parent
                if (!found2) {
                    if (gene_indices2[j] > i) {
                        crossover_inverse_sequence2[i] = crossover_inverse_sequence2[i] + 1;
                    }
                    if (gene_indices2[j] == i) {
                        found2 = true;
                    }
                }
            }
        }
        // crossover afterwards is a simple swap operation with a randomly defined cutoff point
        int cutoff_point = population.GetRandomGeneValue(0);
        for (int i = cutoff_point; i < gene_amount; ++i) {
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
            for (int j = i + 1; j < gene_amount; ++j) {
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
        for (int i = 0; i < gene_amount; ++i) {
            gene_indices1[crossover_position_sequence1[i]] = i;
            gene_indices2[crossover_position_sequence2[i]] = i;
        }

        // transform crossedover continuous value vector back into the explicit city index form
        for (int i = 0; i < gene_amount; ++i) {
            auto it1 = allowed_values.begin();
            std::advance(it1, gene_indices1[i]);
            genes1[i] = *it1;
            auto it2 = allowed_values.begin();
            std::advance(it2, gene_indices2[i]);
            genes2[i] = *it2;
        }
    };
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
            // TODO: do not swap starting point if provided
            // make sure we dont swap in place
            /* what was that? while for one loop, using idx 1 to size?
            while (swapIndex == j) {
                ++swapIndex;
                if (swapIndex == variable_amount + 1) swapIndex = 1;
            }*/
            if (swapIndex == j) {
                ++swapIndex;
                if (swapIndex == variable_amount) swapIndex = 0;
            }
            std::swap(genes[j], genes[swapIndex]);
        }
    }
}

void print_path(const Genotype<VisitedCity> &path,
                const std::map<int, std::vector<const City *>> &node_paths) {
    for (size_t i = 0; i < path.genes.size() - 1; i++) {
        auto from = path.genes[i];
        std::cout << cities[from].name << " ->\n\t\t(";
        auto to = path.genes[i + 1];

        bool reverse = from > to;
        if (reverse) {
            std::swap(from, to);
        } // assure ascending direction to find in distance map
        int key = (from << KEY_UPPER | to);

        if (node_paths.at(key).size() > 0) {
            if (reverse) {
                for (int i = node_paths.at(key).size() - 1; i >= 0; i--) {
                    std::cout << node_paths.at(key)[i]->name << " ";
                }
            } else {
                for (auto c : node_paths.at(key)) {
                    std::cout << c->name << " ";
                }
            }
        }
        std::cout << ")\n";
    }
    std::cout << cities[path.genes[path.genes.size() - 1]].name << std::endl;
}

int main(int argc, char *argv[]) {
    size_t iterations = 500;    // number of generations
    size_t populationSize = 20; // how many different routes are used
    int mutationRate = 5;       // percentage of which a mutation should happen
    int crossoverRate = 70;     // percentage of how often crossover should happen
    bool circle = false;        // flag wheter traveling from last good back to first should happen
    std::set<int> goods;        // provided goods to buy on tour
    int startIdx = -1; // starting point index for route, if not provided just start from first good
    std::string startName = ""; // starting point cityname
    std::string dataPath = "inputs\\romaniaroads.pl";

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

    if (cities.size() < 3) {
        std::cout << "no valid city network found, please contact your local administrator."
                  << std::endl;
        return 0;
    }

    // first check for starting node -> prepend to cluster if provided
    if (false == startName.empty()) {
        if (cityOrdinals.count(startName)) {
            startIdx = cityOrdinals[startName];
        } else {
            std::cout << "given starting point name \"" << startName << "\" invalid!" << std::endl;
            return 0;
        }
    }

    // std::map<int, std::vector<size_t>> cluster;
    std::set<VisitedCity> cluster;
    if (-1 != startIdx) {
        if (cities.size() > (size_t)startIdx) {
            cluster.insert((VisitedCity)startIdx);
        } else {
            std::cout << "starting point index \"" << startIdx << "\" is invalid! "
                      << "range[0 - " << (cities.size() - 1) << "]" << std::endl;
            return 0;
        }
    }

    std::string str = "";
    while (goods.size() < 3 - (-1 != startIdx)) {
        std::cout << "not enough goods provided to justify some work: "
                  << " min: 3, given: " << (goods.size() + (-1 != startIdx)) << std::endl
                  << "please provide indices for cities (0 - " << (cities.size() - 1)
                  << "): <idx1> <idx2> ..." << std::endl;

        std::getline(std::cin, str);
        for (const std::string &idx_city : split(str, " ")) {
            goods.insert(std::stoi(idx_city));
        }
    }

    // select cluster nodes which provide given goods
    for (size_t idx_city = 0; idx_city < cities.size(); ++idx_city) {
        if (std::find(goods.begin(), goods.end(), idx_city) != goods.end()) {
            cluster.insert(idx_city);
        }
    }

    if (cluster.size() < 3) {
        std::cout << "not enough valid cities provided to justify some work: " << cluster.size()
                  << " min: 3" << std::endl;
        return 0;
    }

#ifdef VERBOSE
    std::cout << "selected cluster nodes (" << cluster.size() << "):" << std::endl;
    for (const VisitedCity &node : cluster) {
        std::cout << node << " - " << cities[node].name << std::endl;
    }
#endif

    // calculate distance for cluster nodes (could be done hardcoded because we only use the
    // romaniaroads data file)
    std::map<int, int> node_distances; // key is bitwise combination of node indices, only
                                       // storing one direction (smaller index to bigger)

    std::map<int, std::vector<const City *>>
        node_paths; // key is bitwise combination of node indices, only
                    // storing one direction (smaller index to bigger)

    for (const VisitedCity &from : cluster) {
        for (const VisitedCity &to : cluster) {
            if (from >= to) continue;

            int key = (from << KEY_UPPER | to);
            std::vector<const City *> path;
            int distance = minRoute(from, to, path);
            node_distances.insert(std::make_pair(key, distance));
            node_paths.insert(std::make_pair(key, path));
        }
    }
    std::cout << std::endl << std::endl << std::endl;

#ifdef VERBOSE
    std::cout << "calculated distances (" << node_distances.size()
              << "):" << std::endl; // should be "factorial for sums" - cities with the same good
    for (const auto &dist : node_distances) {
        std::cout << "  " << (dist.first >> KEY_UPPER) << " - " << (dist.first & KEY_LOWER) << " = "
                  << dist.first << ": " << cities[dist.first >> KEY_UPPER].name << " - "
                  << cities[dist.first & KEY_LOWER].name << ": " << dist.second << std::endl;
    }
#endif

    // adjusting size according to variable count
    size_t allowed_values_count = cluster.size();
    crossover_inverse_sequence1.resize(allowed_values_count);
    crossover_inverse_sequence2.resize(allowed_values_count);
    crossover_position_sequence1.resize(allowed_values_count);
    crossover_position_sequence2.resize(allowed_values_count);
    gene_indices1.resize(allowed_values_count);
    gene_indices2.resize(allowed_values_count);

    Population<VisitedCity> population(cluster, iterations, populationSize, crossoverRate / 100.0f,
                                       mutationRate / 100.0f);

    // Fitnesses are negative, so 0 (meaning no distance is traversed at all) would be ideal
    std::cout << "starting evolution with " << iterations << " generations..." << std::endl;
    population.evolve(get_initialization_function(cluster, startIdx),
                      get_evaluation_function(node_distances, circle), mutate_genotype,
                      get_crossover_function(cluster), 0);

    population.print_result();
    auto path = population.GetBestGenotype();
    std::cout << "Path: " << std::endl;
    print_path(path, node_paths);
    return 0;
}
