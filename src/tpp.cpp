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

typedef std::vector<size_t> CitySequence;

void print_sequence(CitySequence &values) {
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

        return sum;
    };
}

void initialize_genotype(CitySequence &genes) {
}

// Declared here so they're re-used in each crossover sequenc
std::vector<int> crossover_inverse_sequence1;
std::vector<int> crossover_inverse_sequence2;
std::vector<int> crossover_position_sequence1;
std::vector<int> crossover_position_sequence2;
void crossover_genotypes(CitySequence &genes1, CitySequence &genes2) {
    // calculate inversion sequences
    int gene_amount = genes1.size();
    // inverse sequence counts the amount of bigger values to the left of the current value and
    // stores it at the corresponding index
    // so for example for the gene [2 3 4 5 1 6 7 8 9] we start with the value 1
    // it is at index 5 and there are 4 bigger values to the left
    // so the value in the inverse sequence at index 1 would be 4
    // this uniquely identifies the gene but also allows to cross with out producing duplicates
    // because translating back a inverse sequence only allows the creation of different values for
    // each gene variable

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
    int cutoff_point = 1; // FIXME: Should be `population.GetRandomGeneValue(0) - 1;` but I want to
                          // avoid a global population
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
        genes1[crossover_position_sequence1[i]] = i + 1;
        genes2[crossover_position_sequence2[i]] = i + 1;
    }
}

void mutate_genotype(CitySequence &genes, double mutation_probability) {
    int variable_amount = genes.size();

    // iterating over each gene
    for (int j = 0; j < variable_amount; ++j) {
        // random chance to mutate gene
        double x = 0.0; // FIXME: `population.GetRandomNormalizedDouble();`

        if (x < mutation_probability) {
            int swapIndex = 1; // FIXME: `population.GetRandomGeneValue(0) - 1;`
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
        // available goods for romaniaroads are (7): 4, 5, 6, 7, 8, 9, 13
    }

    if (goods.size() < 2) {
        std::cout << "not enough goods provided to justify some work: " << goods.size() << " min: 2"
                  << std::endl;
        return 0;
    }

    // 2. read given data containing cities and streets
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    // const std::vector<city> &cities = readInput("data/romaniaroads.pl");
    readInput(dataPath);

    size_t idx = 0;
#ifdef VERBOSE
    std::cout << "cities with connected roads (" << cities.size() << "): " << std::endl;
    for (const auto &city : cities) {
        std::cout << "  " << idx++ << ": " << city.name << " (good: " << city.good
                  << ") with roads (" << city.roads.size() << "): " << std::endl;
        for (Road road : city.roads) {
            std::cout << "\t" << cities[road.to].name << ": " << road.distance << std::endl;
        }
    }
#endif

    // select cluster nodes which provide given goods
    idx = 0;
    std::map<int, std::vector<size_t>> cluster;
    for (const auto &city : cities) {
        if (std::find(goods.begin(), goods.end(), city.good) != goods.end()) {
            if (cluster.count(city.good)) {
                cluster[city.good].push_back(idx);
            } else {
                cluster.insert(std::make_pair(city.good, std::vector<size_t>{idx}));
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
    std::map<int, int> nodeDistances; // key is bitwise combination of node indices, only storing
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
                    nodeDistances.insert(std::make_pair(key, distance));
                }
            }
        }
    }

#ifdef VERBOSE
    std::cout << "calculated distances (" << nodeDistances.size()
              << "):" << std::endl; // should be "factorial for sums" - cities with the same good
    for (const auto &dist : nodeDistances) {
        std::cout << "  " << (dist.first >> KEY_UPPER) << " - " << (dist.first & KEY_LOWER) << " = "
                  << dist.first << ": " << cities[dist.first >> KEY_UPPER].name << " - "
                  << cities[dist.first & KEY_LOWER].name << ": " << dist.second << std::endl;
    }
#endif

    // 3. generate initial purchase order
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    // std::vector<std::unordered_set<size_t>> population; // should use a set for intentional
    // reasoning, but then all the shuffle algorithms would not work :/
    std::vector<std::vector<size_t>> population;
    // TODO: use rng generator class
    static bool init = false;
    if (false == init) {
        init = true;
        srand((size_t)time(NULL));
    }

    for (size_t i = 0; i < populationSize; ++i) {
        // std::unordered_set<size_t> genome;
        std::vector<size_t> genome;
        for (const auto &entry : cluster) {
            // choose random index for multiple options
            int idx = (entry.second.size() == 1) ? 0 : rand() % entry.second.size();
            // genome.insert(entry.second[idx]);
            genome.push_back(entry.second[idx]);
        }

        // population.push_back(shuffle(genome));
        std::random_shuffle(genome.begin() + (-1 != startIdx),
                            genome.end()); // don't shuffle starting point if given
        population.push_back(genome);
    }

    // TEST for all given goods:
    // population.push_back({ 14, 12, 19, 18, 11, 10, 4 }); // 1448
    // population.push_back({ 0, 1, 2, 10, 12, 14, 15 }); // 751
    // population.push_back({ 2, 10, 11, 12, 14, 17, 18 }); // 597
    // populationSize = population.size();

#ifdef VERBOSE
    std::cout << "initial population (" << populationSize << "):" << std::endl;
    for (const auto &pop : population) {
        std::cout << "  ";
        for (int idx : pop) {
            std::cout << idx << ", ";
        }
        std::cout << std::endl;
    }
#endif

    // 4. calculate amount of given iterations (or stop if optimum was found)
    //      - calculate fitness
    //      - select routes with best fitness       -> TODO
    //      - crossover best routes                 -> TODO
    //      - mutate a small percentage             -> TODO
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    std::vector<int> fitness(populationSize); // TODO: should be property of population entry
    std::cout << "TODO: handle mutation rate " << mutationRate << std::endl;
    std::cout << "iterating " << iterations << " times..." << std::endl;
    int min = INT16_MAX;
    std::vector<size_t> superiorAlphaGenom;
    for (size_t i = 0; i < iterations; ++i) {
        for (size_t p = 0; p < populationSize; ++p) {
            // calculating fitness
            int sum = 0, from, to;
            std::vector<size_t>::iterator it_to;
            for (auto it_from = population[p].begin(); it_from != population[p].end(); ++it_from) {
                from = *it_from;
                it_to = std::next(it_from, 1);
                if (it_to == population[p].end()) {
                    if (false == circle) break;

                    // if purchaser has to travel back to start add closest distance
                    it_to = population[p].begin();
                }

                to = *it_to;
                if (from > to)
                    std::swap(from, to); // assure ascending direction to find in distance map

                sum += nodeDistances[(from << KEY_UPPER | to)];
            }
            fitness[p] = sum;
            if (sum < min) {
                min = sum;
                superiorAlphaGenom = population[p];
            }
        }
    }

    // 5. output best route
    // --------------------------------------------------------------------------------------------------------------------------------------------------------
    /*std::cout << "calculated fitnesses:" << std::endl;
    min = fitness[0];
    for (int fit : fitness)
    {
        if (fit < min) min = fit;
        std::cout << fit << std::endl;
    }*/
    std::cout << "minimum distance: " << min << std::endl;
    std::cout << "for genome: ";
    for (int chromosome : superiorAlphaGenom) {
        std::cout << chromosome << ", ";
    }
    if (circle) {
        std::cout << superiorAlphaGenom[0];
    }
    std::cout << std::endl;
    return 0;
}