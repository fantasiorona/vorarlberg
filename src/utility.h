#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
//#include <sstream>
#include <vector>

struct Road;

// marketplace
struct City {
    std::string name;
    int good;                // only one good -> calculated from name length
    std::vector<Road> roads; // TODO: could use pointers for symetrical connections

    City() {
    }

    City(std::string name) : name(name), good(name.length()) {
    }
};

// symmetrical connection between cities
struct Road {
    int to; // index from cities vector
    int distance;

    Road(int to, int distance) : to(to), distance(distance) {
    }
};

std::map<std::string, int> cityOrdinals; // mapping name to vector index
std::vector<City> cities;

// TODO: add actual visited cities to path
int minRoute(const City *a, const City *b, std::vector<const City *> &path) {
    int minRoute = 0;                  // a->good;
    std::vector<const City *> visited; // mark already visited nodes
    std::multimap<int, const City *>
        weights; // store distances in sorted map, allows duplicate keys
    do {
        visited.push_back(a);
        for (const auto &road : a->roads) {
            if (std::find(visited.begin(), visited.end(), &cities[road.to]) == visited.end()) {
                // put min distances on heap
                // weights.insert(std::make_pair(minRoute + road.distance, &cities[road.to])); //
                // without buying goods on visited cities
                weights.insert(std::make_pair(minRoute + road.distance + cities[road.to].good,
                                              &cities[road.to]));
            }
        }

        do {
            // remove from heap, update minRoute
            minRoute = weights.begin()->first;
            a = weights.begin()->second;
            // weights.erase(minRoute); // erases all keys
            weights.erase(weights.find(minRoute));
        } while (std::find(visited.begin(), visited.end(), a) != visited.end());
    } while (a != b);

    return minRoute;
}

std::vector<std::string> split(const std::string &s, const std::string &delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::vector<std::string> res;
    std::string token;
    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

const std::vector<City> &readInput(const std::string &path) {
    if (!cities.empty()) cities.clear();

    std::cout << "reading data from file " << path << std::endl;
    std::ifstream file(path);
    std::string line;
    std::string str_from, str_to;
    while (getline(file, line)) {
        if (line.find('%') == 0) { // skipping comments
            continue;
        } else if (line.length() < 3) { // done at first empty line, only contains newline
            break;
        }

        // parsing "road(arad, zerind, 75)."
        auto start = line.find('(');
        auto end = line.find(')');
        std::string connection = line.substr(start + 1, end - start - 1);
        // splitting "arad, zerind, 75"
        std::vector<std::string> data = split(connection, ", ");
        if (data.size() != 3) {
            std::cout << "could not parse line: " << line << std::endl;
            throw std::invalid_argument("invalid data provided");
        }

        str_from = data[0];
        str_to = data[1];
        // although inserting an already existing key into a map does nothing
        // it is not needed to construct a new object
        if (cityOrdinals.count(str_from) == 0) {
            cities.push_back(City(str_from));
            cityOrdinals.insert(std::make_pair(str_from, cities.size() - 1));
        }

        if (cityOrdinals.count(str_to) == 0) {
            cities.push_back(City(str_to));
            cityOrdinals.insert(std::make_pair(str_to, cities.size() - 1));
        }

        Road roadTo(cityOrdinals[str_to],
                    std::stoi(data[2])); // stoi throws an invalid_argument exception if no
                                         // conversion could be performed
        Road roadFrom(cityOrdinals[str_from], std::stoi(data[2]));

        cities[cityOrdinals[str_from]].roads.push_back(roadTo);
        cities[cityOrdinals[str_to]].roads.push_back(roadFrom);
    }
    return cities;
}