#ifndef CPU_PREDEFINED_TYPES_H
#define CPU_PREDEFINED_TYPES_H

/**
 * Data structure that contains pre calced product for pair of matrices i.e
 * For matrices p,q of size n it contains its product
 */
typedef std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>> PrecalcMap;


#endif //CPU_PREDEFINED_TYPES_H
