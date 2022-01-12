
#ifndef CPU_PARSERS_H
#define CPU_PARSERS_H


#include <cmath>
#include <algorithm>
#include <unordered_set>


#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "semi_local.h"


int *split(std::string str, const std::string &token, int arr_length) {
    int *result = new int[arr_length];
    int i = 0;
    while (str.size()) {
        int index = str.find(token);
        if (index != std::string::npos) {
            result[i] = stoi(str.substr(0, index));
            i++;
            str = str.substr(index + token.size());
        } else {
            str = "";
        }
    }
    return result;
}


std::pair<int, std::pair<std::string, std::string>> parse_input_file(const std::string &filename) {
    std::ifstream file(filename);
    std::string name;
    std::string sequence;
    std::string length;
    std::getline(file, name);
    std::getline(file, length);
    std::getline(file, sequence);
    return std::make_pair(stoi(length), std::make_pair(name, sequence));
};





#endif //CPU_PARSERS_H
