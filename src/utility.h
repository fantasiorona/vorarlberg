#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
//#include <sstream>
#include <vector>

struct road;

struct city // marketplace
{
    std::string name;
    int good; // only one good -> calculated from name length
    std::vector<road> roads; // TODO: could use pointers for symetrical connections

    city ()
    {
    }

    city (std::string name) : name(name), good(name.length())
    {
    }
};

// roads are symetrical
struct road // connection between cities
{
    int to; // index from cities vector
    int distance;

    road (int to, int distance) : to(to), distance(distance)
    {
    }
};

std::map<std::string, int> cityOrdinals; // mapping name to vector index
std::vector<city> cities;

int minRoute (const city *a, const city *b)
{
    int minRoute = 0;
    std::vector<const city *> visited; // mark already visited nodes
    std::multimap<int, const city* > weights; // store distances in sorted map, allows duplicate keys
    do
    {
        visited.push_back(a);
        for (const auto &road : a->roads)
        {
            if (std::find(visited.begin(), visited.end(), &cities[road.to]) == visited.end())
            {
                weights.insert(std::make_pair(road.distance + minRoute, &cities[road.to])); // put on heap
            }
        }

        do
        {
            // remove from heap, update minRoute
            minRoute = weights.begin()->first;
            a = weights.begin()->second;

            //weights.erase(minRoute); // erases all keys
            weights.erase(weights.find(minRoute));
        }
        while (std::find(visited.begin(), visited.end(), a) != visited.end());
    }
    while (a != b);

    return minRoute;
}

std::vector<std::string> split (const std::string &s, const std::string &delimiter)
{
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::vector<std::string> res;
    std::string token;
    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos)
    {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

const std::vector<city> & readInput (const std::string &path)
{
    if (!cities.empty()) cities.clear();

    std::cout << "reading data from file " << path << std::endl;
    std::ifstream file(path);
    std::string line;
    std::string str_from, str_to;
    while (getline(file, line))
    {
        std::cout << "parsing line " << line;
        if (line.find('%') == 0) // skipping comments
        {
            continue;
        }
        else if (line.length() < 3) // done at first empty line, only contains newline
        {
            break;
        }

        // parsing "road(arad, zerind, 75)."
        auto start = line.find('(');
        auto end = line.find(')');
        std::string connection = line.substr(start + 1, end - start - 1);
        // splitting "arad, zerind, 75"
        std::vector<std::string> data = split(connection, ", ");
        if (data.size() != 3)
        {
            std::cout << "could not parse line: " << line << std::endl;
            throw std::invalid_argument("invalid data provided");
        }

        str_from = data[0];
        str_to = data[1];
        // although inserting an already existing key into a map does nothing
        // it is not needed to construct a new object
        if (cityOrdinals.count(str_from) == 0)
        {
            cities.push_back(city(str_from));
            cityOrdinals.insert(std::make_pair(str_from, cities.size() - 1));
        }

        if (cityOrdinals.count(str_to) == 0)
        {
            cities.push_back(city(str_to));
            cityOrdinals.insert(std::make_pair(str_to, cities.size() - 1));
        }

        road roadTo(cityOrdinals[str_to], std::stoi(data[2])); // stoi throws an invalid_argument exception if no conversion could be performed
        road roadFrom(cityOrdinals[str_from], std::stoi(data[2]));

        cities[cityOrdinals[str_from]].roads.push_back(roadTo);
        cities[cityOrdinals[str_to]].roads.push_back(roadFrom);
    }
    return cities;
}