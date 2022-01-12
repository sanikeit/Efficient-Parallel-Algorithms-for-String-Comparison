#ifndef CPU_PERMUTATIONS_ENCODING_H
#define CPU_PERMUTATIONS_ENCODING_H

#include <cmath>
#include "matrices.h"

/**
 * We encode permutation matrices naively. It only works till specified size (around 7-8).
 * Idea is simple, we embedded cols position in each row of non-zero entires in 32(64) bit word
 */
namespace std {

    template<>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const std::pair<int, int> &pair) const {
            return hash<long long>()(((long long) pair.first) ^ (((long long) pair.second) << 32));
        }
    };


    template<>
    struct hash<AbstractPermutation> {
        std::size_t operator()(const AbstractPermutation &k) const {
            using std::size_t;
            using std::hash;
            using std::string;

            size_t sum = 0;
            auto bits_per_symbol = int(std::ceil(log2(k.row_size)));

            for (int i = 0; i < k.row_size; ++i) {
                sum = (sum << bits_per_symbol) + k.get_col_by_row(i);

            }
            return sum;
        }

    };

    bool operator==(const Permutation &p1, const Permutation &p2) {
        return std::hash<AbstractPermutation>()(p1) == std::hash<AbstractPermutation>()(p1);
    }

}


#endif //CPU_PERMUTATIONS_ENCODING_H
