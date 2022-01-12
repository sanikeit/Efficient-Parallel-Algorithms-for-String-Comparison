#ifndef CPU_SEQUENCE_GENERATORS_H
#define CPU_SEQUENCE_GENERATORS_H


#include <vector>
#include <cstdlib>

/**
 * Generate sequence of int numbers
 * @param size int
 * @param alphabet_size int, default = 4
 * @return vector<int> of size @size
 */
template<class T>
std::vector<T> gen_vector_seq(int size, int alphabet_size = 4) {
    auto v = std::vector<T>();
    v.reserve(size);
    for (int i = 0; i < size; ++i) {
        v.push_back(T(rand() % alphabet_size));
    }
    return v;
}

/**
 * Generate sequence of int numbers
 * @param size int
 * @param alphabet_size int, default = 4
 * @return vector<int> of size @size
 */
int * gen_array(int size, int alphabet_size = 4) {
    int * a = new int[size];
    for (int i = 0; i < size; ++i) {
        a[i] = (rand() % alphabet_size);
    }
    return a;
}



#endif //CPU_SEQUENCE_GENERATORS_H
