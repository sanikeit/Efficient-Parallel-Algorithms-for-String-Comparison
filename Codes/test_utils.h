#ifndef CPU_TEST_UTILS_H
#define CPU_TEST_UTILS_H

const long long R = 4294967279;
const long long M = 4294967291;

//static const int length = 1024*8;
long long hash(AbstractPermutation &arr, int size) {
    long long hash = 0;
    for (int i = 0; i < size; i++) {
        hash = (R * hash + arr.get_row_by_col(i)) % M;
    }
    return hash;
}



#endif //CPU_TEST_UTILS_H
