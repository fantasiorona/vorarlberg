#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
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

// modifyed djikstra that also stores the path information
// based on:
// https://www.codingame.com/playgrounds/1608/shortest-paths-with-dijkstras-algorithm/dijkstras-algorithm
int minRoute(size_t a, size_t b, std::vector<const City *> &path) {
    // keeps track of already visited cities
    std::vector<bool> visited(cities.size(), false);
    // keeps track of current minimum distances to all cities
    std::vector<int> distances(cities.size(), 9999999);
    // keeps track of current minimum paths
    std::vector<std::vector<const City *>> paths(cities.size());
    // keeps track of cities we have to visit
    std::vector<size_t> to_visit;

    // starting at a
    to_visit.push_back(a);
    // distance to start city is 0
    distances[a] = 0;
    size_t current_city = a;
    // do the loop while we have no more cities to visit or we have reached b
    while (!to_visit.empty() && current_city != b) {
        // distance to next nearest city we can visit, set basically to infinity at beginning
        int nearest_city_distance = 9999999;
        size_t next_city = 0;
        // iterate over all city to visity
        for (size_t i = 0; i < to_visit.size(); i++) {
            // comparing with currently smalles distance to find the next city
            if (distances[to_visit[i]] < nearest_city_distance) {
                next_city = i;
                nearest_city_distance = distances[to_visit[i]];
            }
        }
        current_city = to_visit[next_city];
        // remove current_city from cities to visit
        std::swap(to_visit[next_city], to_visit[to_visit.size() - 1]);
        to_visit.resize(to_visit.size() - 1);

        // iterate over all neighbors of current city
        for (auto road : cities[current_city].roads) {
            // and calculate the total distance to the start city through the current path
            int distance = distances[current_city] + road.distance + cities[road.to].good;
            // if we found a shorter path
            if (distances[road.to] > distance) {
                // store new shortest distance
                distances[road.to] = distance;
                // store the path:
                // clear previous path
                paths[road.to].clear();
                // store path to previous city
                for (auto c : paths[current_city]) {
                    paths[road.to].push_back(c);
                }
                // add current city
                paths[road.to].push_back(&cities[road.to]);
            }
            // store cities we have not yet visited
            if (!visited[road.to]) {
                to_visit.push_back(road.to);
            }
        }
        visited[current_city] = true;
    }

    // set return path if available
    if (paths[b].size() > 0) {
        for (size_t i = 0; i < paths[b].size() - 1; i++) {
            path.push_back(paths[b][i]);
        }
    }
    return distances[b];
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