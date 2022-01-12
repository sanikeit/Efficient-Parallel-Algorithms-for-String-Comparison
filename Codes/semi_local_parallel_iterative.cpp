#include <string>
#include <iostream>
#include <chrono>
#include "semi_local.h"
#include "parsers.h"
#include "sequence_generators.h"
#include "test_utils.h"

#include "distribution_matrice.h"

int main(int argc, char *argv[]) {
    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    int a_size = name_content_a.first;
    int b_size = name_content_b.first;
    auto a_name = name_content_a.second.first;
    auto b_name = name_content_b.second.first;
    int * a = split(name_content_a.second.second,",",a_size);
    int * b = split(name_content_b.second.second,",",b_size);

    auto perm = Permutation(a_size+b_size,a_size+b_size);

    auto beg = std::chrono::high_resolution_clock::now();
    semi_local::sticky_braid_mpi<int, false, true>(perm, a, a_size, b, b_size, thds);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
    std::cout << 0  << std::endl; // some preprocess
    std::cout << "elapsed_time  " << elapsed_time  << std::endl; // algo time
    std::cout << hash(perm, perm.row_size) << std::endl;

    
    // std::vector<int> perm_(a_size + b_size);
    // for(int i = 0; i < perm.row_size; ++i){
    //     perm_[i] = perm.get_row_by_col(i);
    //     std::cout << perm_[i] << " ";
    // }
    // auto H = semi_local_lcs(perm_,a_size);
    


    std::cout<< a_size<<std::endl;
    std::cout<< b_size<<std::endl;
    std::cout<< a_name<<std::endl;
    std::cout<< b_name;

    delete[] a;
    delete[] b;
}