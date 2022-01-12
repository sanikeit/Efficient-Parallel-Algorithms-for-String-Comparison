#ifndef CPU_STEADY_ANT_H
#define CPU_STEADY_ANT_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <climits>
#include <omp.h>
#include "matrices.h"
#include "permutations_encoding.h"
#include "dominance_sum_queries.h"
#include "semi_local.h"
#include "predefined_types.h"


namespace distance_unit_monge_product {

    /**
     * Contains naïve implementation of unit Monge distance product.
     * First, for each matrix a dominance sum matrices are constructed.
     * Then O(n^3) time complexity  distance product of explicit matrices is performed
     */
    namespace naive {
        bool top_right_summator(int cur_row, int cur_col, int row_bound, int col_bound) {
            return cur_row >= row_bound && cur_col < col_bound;
        }


        void mult_dist(AbstractPermutation *m, AbstractPermutation *n, AbstractPermutation *distance_product) {
            auto dominance_m = Matrix(m->row_size + 1, m->col_size + 1);
            get_dominance_matrix(*m, top_right_summator, &dominance_m);
            auto dominance_n = Matrix(n->row_size + 1, n->col_size + 1);
            get_dominance_matrix(*n, top_right_summator, &dominance_n);

            auto row_size = (m->row_size + 1);
            auto col_size = (n->col_size + 1);
            auto dominance_c = Matrix(row_size, col_size);

            for (int i = 0; i < row_size; ++i) {
                for (int k = 0; k < col_size; ++k) {
                    auto tmp = INT_MAX;
                    for (int j = 0; j < (m->col_size + 1); ++j) {
                        tmp = std::min(dominance_m.get_element_at(i, j) + dominance_n.get_element_at(j, k), tmp);
                    }
                    dominance_c.set_element_at(i, k, tmp);
                }
            }

            for (int i = 0; i < m->row_size; ++i) {
                for (int j = 0; j < n->col_size; ++j) {
                    auto cross_diff =
                            (dominance_c.get_element_at(i, j + 1) + dominance_c.get_element_at(i + 1, j)) -
                            (dominance_c.get_element_at(i, j) + dominance_c.get_element_at(i + 1, j + 1));

                    if (cross_diff == 1) distance_product->set_point(i, j);

                }
            }

        }
    };



    /**
     * Contains sophisticated implementations of unit Monge distance product.
     * It uses recent successes over unit monge distance multiplication.
     * Namespace contains contains implementation of so called steady ant algorithm that performs O(nlogn) time complexity
     * multiplication of two permutation matrices.
     */
    namespace steady_ant {

        namespace _details {

            /**
            * see theorem 5.21
            * Allows get P_{b,a} when you have P_{a,b}
            */
            void fill_permutation_ba(AbstractPermutation *ab, AbstractPermutation *ba, int m, int n) {
                ba->unset_all();
                for (int i = 0; i < ab->row_size; ++i) {
                    auto col = ab->get_col_by_row(i);
                    if (col != NOPOINT) ba->set_point(n + m - 1 - i, m + n - 1 - col);
                }
            }


            /**
             * Given squared permutation matrix p_{i} get slice p[,start_inclusive:end_exclusive] and map it to new coordinates
             * to get new squared matrix of size end_exclusive-start_inclusive and mapping of row coordinates (new to old)
             * @param p_i squared permutation matrix
             * @param col_start_inclusive inclusive index of start of slice
             * @param col_end_exclusive exclusive index of end of slice
             * @return mapping of row_{i+1} -> row_{i} and p_{i+1}, mapping of col_{i+1} -> col_{i} implicit via offset start_inclusive
             */
            void get_vertical_slice(AbstractPermutation *p_i, int col_start_inclusive, int col_end_exclusive,
                                    int *cur_row_to_prev_row_mapping, AbstractPermutation *succ_pi) {

                auto ptr = 0;
                for (int row = 0; row < p_i->row_size; ++row) {
                    auto col = p_i->get_col_by_row(row);
                    if (col >= col_start_inclusive && col < col_end_exclusive) {
                        cur_row_to_prev_row_mapping[ptr] = row;
                        succ_pi->set_point(ptr, col - col_start_inclusive);
                        ptr++;
                    }
                }
            }


            /**
             * Given squared permutation matrix q_{i} get slice p[start_inclusive:end_exclusive,:] and map it to new coordinates
            * to get new squared matrix of size end_exclusive-start_inclusive and mapping of col coordinates (new to old)
            * @param q_i squared permutation matrix
            * @param row_start_inclusive inclusive index of start of slice
            * @param row_end_exclusive exclusive index of end of slice
            * @return mapping of col_{i+1} -> col_{i} and p_{i+1}, mapping of row_{i+1} -> row_{i} implicit via offset start_inclusive
            */
            void get_horizontal_slice(AbstractPermutation *q_i, int row_start_inclusive, int row_end_exclusive,
                                      int *cur_col_to_prev_col_mapping, AbstractPermutation *succ_pi) {

                auto ptr = 0;
                for (int col = 0; col < q_i->col_size; ++col) {
                    auto row = q_i->get_row_by_col(col);
                    if (row >= row_start_inclusive && row < row_end_exclusive) {
                        cur_col_to_prev_col_mapping[ptr] = col;
                        succ_pi->set_point(row - row_start_inclusive, ptr);

                        ptr++;
                    }
                }
            }

            /**
             * Maps non-zero entries of shrinked n/2Xn/2 permutation matrix to nXn matrix. Aka maps positions of non-zero
             * elements of matrix shrinked to flattened matrix  according to row_mapper and col_mapper
             * @param shrinked
             * @param row_mapper
             * @param col_mapper
             * @param flattened
             */
            inline void
            inverse_mapping(AbstractPermutation *shrinked, const int *row_mapper, const int *col_mapper,
                            AbstractPermutation *flattened) {
                //could be parallelized
                flattened->unset_all();

                //#pragma omp parallel for
                for (int cur_col = 0; cur_col < shrinked->col_size; ++cur_col) {
                    auto old_col = col_mapper[cur_col]; /// consecutive access

                    auto cur_row = shrinked->get_row_by_col(cur_col); // consecutive access
                    auto old_row = row_mapper[cur_row]; // non consecutive access
                    flattened->set_point(old_row, old_col);
                }
            }

            /**
             * Maps non-zero entries of shrinked n/2Xn/2 permutation matrix to nXn matrix. Aka maps positions of non-zero
             * elements of matrix shrinked to flattened matrix  according to row_mapper and col_mapper
             * @param shrinked
             * @param row_mapper
             * @param col_mapper
             * @param flattened
             */
            inline void
            inverse_mapping_without_pre_clearing(AbstractPermutation *shrinked, const int *row_mapper,
                                                 const int *col_mapper,
                                                 AbstractPermutation *flattened) {

                //#pragma omp parallel for
                for (int cur_col = 0; cur_col < shrinked->col_size; ++cur_col) {
                    auto old_col = col_mapper[cur_col]; /// consecutive access

                    auto cur_row = shrinked->get_row_by_col(cur_col); // consecutive access
                    auto old_row = row_mapper[cur_row]; // non consecutive access
                    flattened->set_point(old_row, old_col);
                }
            }


            inline void
            ant_passage(AbstractPermutation *r_lo, AbstractPermutation *r_hi, int n, std::vector<int> &good_row_pos,
                        std::vector<int> &good_col_pos) {

                //ant passage
                auto end_row = -1;
                auto end_col = r_lo->col_size + 1;
                auto cur_row = r_hi->row_size;
                auto cur_col = -1;

                auto rhi = 0;
                auto rlo = 0;

                bool is_went_right = false; // went up
                while (true) {

                    if (end_col == cur_col && end_row == cur_row) break;

                    if (cur_row == 0) break;
                    //TODO is new
                    if (cur_col == n) break;
                    //
                    auto dominance_row = cur_row - 1;
                    auto dominance_col = cur_col + 1;

                    //prev step
                    if (is_went_right) {
                        rhi = dominance_sum_counting::bottom_right_arrow::right_move(dominance_row,
                                                                                     dominance_col - 1, rhi,
                                                                                     r_hi);
                        rlo = dominance_sum_counting::top_left_arrow::right_move(dominance_row,
                                                                                 dominance_col - 1, rlo,
                                                                                 r_lo);
                    } else {
                        rhi = dominance_sum_counting::bottom_right_arrow::up_move(dominance_row + 1,
                                                                                  dominance_col, rhi,
                                                                                  r_hi);
                        rlo = dominance_sum_counting::top_left_arrow::up_move(dominance_row + 1,
                                                                              dominance_col,
                                                                              rlo, r_lo);
                    }

                    if (rhi - rlo < 0) {
                        is_went_right = true;
                        cur_col++;
                    } else if (rhi - rlo == 0) {
                        is_went_right = false;
                        cur_row--;
                    } else {
                        std::cout << "Impissble" << std::endl;
                    }

                    if (dominance_col > 0) {
                        auto delta_above_left =
                                dominance_sum_counting::bottom_right_arrow::left_move(dominance_row,
                                                                                      dominance_col, rhi,
                                                                                      r_hi) -
                                dominance_sum_counting::top_left_arrow::left_move(dominance_row,
                                                                                  dominance_col,
                                                                                  rlo, r_lo);
                        auto delta_below_right =
                                dominance_sum_counting::bottom_right_arrow::down_move(dominance_row,
                                                                                      dominance_col, rhi,
                                                                                      r_hi) -
                                dominance_sum_counting::top_left_arrow::down_move(dominance_row,
                                                                                  dominance_col,
                                                                                  rlo, r_lo);

                        if (delta_above_left < 0 && delta_below_right > 0) {
                            good_row_pos.push_back(dominance_row);
                            good_col_pos.push_back(dominance_col - 1);

                        }
                    }


                }
                // end ant passage

            }

        }

        /**
         * Non optimized version of steady ant algorithm that allocated new memory on each recursion step and
         * uses  precomputed if had
         * @param p
         * @param q
         * @param map
         * @return
         */
        AbstractPermutation *steady_ant(AbstractPermutation *p, AbstractPermutation *q, PrecalcMap &map) {
            using namespace _details;

            auto n = p->col_size;


            if (n <= map.size()) {
                auto precalced = new Permutation(n, n, map[n][std::hash<AbstractPermutation>()(
                        *p)][std::hash<AbstractPermutation>()(*q)]);
                return precalced;
            }

            if (n == 1) {
                //base case when common dimension is one
                auto row = p->get_row_by_col(0);
                auto col = q->get_col_by_row(0);
                auto matrix = new Permutation(p->row_size, q->col_size);

                delete p;
                delete q;

                if (row != NOPOINT && col != NOPOINT) matrix->set_point(row, col);
                return matrix;
            }


            int spliter = p->col_size / 2;

            auto p_lo_row_mapper = new int[spliter];
            auto p_lo = new Permutation(spliter, spliter);

            get_vertical_slice(p, 0, spliter, p_lo_row_mapper, p_lo);

            auto p_hi_row_mapper = new int[p->col_size - spliter];
            auto p_hi = new Permutation(p->col_size - spliter, p->col_size - spliter);

            get_vertical_slice(p, spliter, p->col_size, p_hi_row_mapper, p_hi);


            auto q_lo_col_mapper = new int[spliter];
            auto q_lo = new Permutation(spliter, spliter);
            get_horizontal_slice(q, 0, spliter, q_lo_col_mapper, q_lo);

            auto r_lo = p; // reuse since we no need to use p further


            auto product = steady_ant(p_lo, q_lo, map);


            inverse_mapping(product, p_lo_row_mapper, q_lo_col_mapper, r_lo);
            delete product;

            auto q_hi_col_mapper = new int[p->col_size - spliter];
            auto q_hi = new Permutation(p->col_size - spliter, p->col_size - spliter);
            get_horizontal_slice(q, spliter, q->row_size, q_hi_col_mapper, q_hi);
            auto r_hi = q; // reuse since we no need to use p further
            product = steady_ant(p_hi, q_hi, map);


            inverse_mapping(product, p_hi_row_mapper, q_hi_col_mapper, r_hi);

            delete product;
            delete[] p_lo_row_mapper;
            delete[] q_lo_col_mapper;
            delete[] p_hi_row_mapper;
            delete[] q_hi_col_mapper;

            std::vector<int> good_row_pos;
            std::vector<int> good_col_pos;
            ant_passage(r_lo, r_hi, n, good_row_pos, good_col_pos);

            // merge r_lo to r_hi
            for (int i = 0; i < n; ++i) {
                auto col = r_lo->get_col_by_row(i);
                if (col == NOPOINT) continue;
                r_hi->set_point(i, col);
            }

            // add good points
            for (int i = 0; i < good_col_pos.size(); ++i) {
                auto col = good_col_pos[i];
                auto row = good_row_pos[i];
                r_hi->set_point(row, col);
            }

            delete r_lo;

            return r_hi;
        }


        /**
         * Optimized version. Uses linear memory. You need to allocate in total 4n for memory block matrix;
         * p uses 2n memory, q uses 2n memory.
         * Result can be obtained from p matrix
         * @param p permutation of size nXn
         * @param q permutation of size nXn
         * @param memory_block_matrices  4n memory of p and q
         * @param free_space_matrices  extra 4n mamory that algorithm uses to split matrices on sub parts
         * @param memory_block_indices  8n memory for stroring memory indices in sequential case
         * @param map precalc for  size at least 1, optimal size is 5
         */
        void steady_ant_optimized_seq(
                AbstractPermutation *p, AbstractPermutation *q, int *memory_block_matrices, int *free_space_matrices,
                int *memory_block_indices, PrecalcMap &map) {
            using namespace _details;

            auto n = p->row_size;

            if (n <= map.size()) {
                auto precalced = PermutationPreAllocated(n, n, memory_block_matrices, memory_block_matrices + n,
                                                         map[n][std::hash<AbstractPermutation>()(
                                                                 *p)][std::hash<AbstractPermutation>()(*q)]);
                return;
            }

            int spliter = n / 2;


            auto p_lo_row_mapper = memory_block_indices;
            auto p_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices, free_space_matrices + spliter);
            get_vertical_slice(p, 0, spliter, p_lo_row_mapper, &p_lo);


            auto q_lo_col_mapper = memory_block_indices + spliter;
            auto q_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices + 2 * spliter,
                                                free_space_matrices + 3 * spliter);
            get_horizontal_slice(q, 0, spliter, q_lo_col_mapper, &q_lo);


            auto p_hi_row_mapper = memory_block_indices + 2 * spliter;
            auto p_hi = PermutationPreAllocated(n - spliter, n - spliter, free_space_matrices + 4 * spliter,
                                                free_space_matrices + 4 * spliter + (n - spliter));
            get_vertical_slice(p, spliter, n, p_hi_row_mapper, &p_hi);


            auto q_hi_col_mapper = memory_block_indices + 2 * spliter + (n - spliter);
            auto q_hi = PermutationPreAllocated(n - spliter, n - spliter,
                                                free_space_matrices + 4 * spliter + 2 * (n - spliter),
                                                free_space_matrices + 4 * spliter + 3 * (n - spliter));
            get_horizontal_slice(q, spliter, q->row_size, q_hi_col_mapper, &q_hi);

            // hack
            auto r_lo = p;
            auto r_hi = q;


            steady_ant_optimized_seq(&p_lo, &q_lo, free_space_matrices, memory_block_matrices,
                                     memory_block_indices + 2 * n, map);


            steady_ant_optimized_seq(&p_hi, &q_hi, free_space_matrices + 4 * spliter,
                                     memory_block_matrices + 4 * spliter,
                                     memory_block_indices + 2 * n, map);

            inverse_mapping(&p_lo, p_lo_row_mapper, q_lo_col_mapper, r_lo);
            inverse_mapping(&p_hi, p_hi_row_mapper, q_hi_col_mapper, r_hi);

            std::vector<int> good_row_pos;
            std::vector<int> good_col_pos;
            ant_passage(r_lo, r_hi, n, good_row_pos, good_col_pos);

            // merge r_hi to r_lo
            for (int i = 0; i < n; ++i) {
                auto col = r_hi->get_col_by_row(i);
                if (col == NOPOINT) continue;
                r_lo->set_point(i, col);
            }

            // add good points
            for (int i = 0; i < good_col_pos.size(); ++i) {
                auto col = good_col_pos[i];
                auto row = good_row_pos[i];
                r_lo->set_point(row, col);
            }


        }


        /*
        * Optimized version of steady ant algorithm that that uses preallocated memory blocks for mapping and
        * matrices itself. Also it uses precomputed values via PrecalcMap.
        * Exactly 8n memory required for stroring matrices.
        * 2nlog(2n) memory is required for parallel version for memory mapping.
        */
        void steady_ant_with_precalc_and_memory(
                AbstractPermutation *p, AbstractPermutation *q, int *memory_block_matrices, int *free_space_matrices,
                int *memory_block_indices, PrecalcMap &map, int total_memory, int nested_parall_regions = 2) {
            using namespace _details;

            auto n = p->row_size;

            if (n <= map.size()) {
                auto precalced = PermutationPreAllocated(n, n, memory_block_matrices, memory_block_matrices + n,
                                                         map[n][std::hash<AbstractPermutation>()(
                                                                 *p)][std::hash<AbstractPermutation>()(*q)]);
                return;
            }

            int spliter = n / 2;


            auto p_lo_row_mapper = memory_block_indices;
            auto p_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices, free_space_matrices + spliter);
            get_vertical_slice(p, 0, spliter, p_lo_row_mapper, &p_lo);


            auto q_lo_col_mapper = memory_block_indices + spliter;
            auto q_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices + 2 * spliter,
                                                free_space_matrices + 3 * spliter);
            get_horizontal_slice(q, 0, spliter, q_lo_col_mapper, &q_lo);


            auto p_hi_row_mapper = memory_block_indices + 2 * spliter;
            auto p_hi = PermutationPreAllocated(n - spliter, n - spliter, free_space_matrices + 4 * spliter,
                                                free_space_matrices + 4 * spliter + (n - spliter));
            get_vertical_slice(p, spliter, n, p_hi_row_mapper, &p_hi);


            auto q_hi_col_mapper = memory_block_indices + 2 * spliter + (n - spliter);
            auto q_hi = PermutationPreAllocated(n - spliter, n - spliter,
                                                free_space_matrices + 4 * spliter + 2 * (n - spliter),
                                                free_space_matrices + 4 * spliter + 3 * (n - spliter));
            get_horizontal_slice(q, spliter, q->row_size, q_hi_col_mapper, &q_hi);

            // now we have small matrices in free space, and p,q may be overwritten
//        free_space_mappings + 2*n

            // hack
            auto r_lo = p;
            auto r_hi = q;

            int on_parts = (total_memory - 2 * n) / 2;


            if (nested_parall_regions > 0) {

#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {

#pragma omp task depend(in: on_parts)
                        steady_ant_with_precalc_and_memory(&p_lo, &q_lo, free_space_matrices, memory_block_matrices,
                                                           memory_block_indices + 2 * n, map, on_parts,
                                                           nested_parall_regions - 1);
#pragma omp task depend(in: on_parts)
                        steady_ant_with_precalc_and_memory(&p_hi, &q_hi, free_space_matrices + 4 * spliter,
                                                           memory_block_matrices + 4 * spliter,
                                                           memory_block_indices + 2 * n + on_parts, map, on_parts,
                                                           nested_parall_regions - 1);

#pragma omp task depend(out: on_parts)
                        {

                            inverse_mapping(&p_lo, p_lo_row_mapper, q_lo_col_mapper, r_lo);
                            inverse_mapping(&p_hi, p_hi_row_mapper, q_hi_col_mapper, r_hi);

                            std::vector<int> good_row_pos;
                            std::vector<int> good_col_pos;
                            ant_passage(r_lo, r_hi, n, good_row_pos, good_col_pos);


                            // merge r_hi to r_lo
                            for (int i = 0; i < n; ++i) {
                                auto col = r_hi->get_col_by_row(i);
                                if (col == NOPOINT) continue;
                                r_lo->set_point(i, col);
                            }

                            // add good points
                            for (int i = 0; i < good_col_pos.size(); ++i) {
                                auto col = good_col_pos[i];
                                auto row = good_row_pos[i];
                                r_lo->set_point(row, col);
                            }


                        }
//#pragma omp taskwait
                    };



                    // new matrix in r_lo

                }

            } else {
                steady_ant_with_precalc_and_memory(&p_lo, &q_lo, free_space_matrices, memory_block_matrices,
                                                   memory_block_indices + 2 * n, map, on_parts, nested_parall_regions);
                steady_ant_with_precalc_and_memory(&p_hi, &q_hi, free_space_matrices + 4 * spliter,
                                                   memory_block_matrices + 4 * spliter,
                                                   memory_block_indices + 2 * n + on_parts, map, on_parts,
                                                   nested_parall_regions);

                inverse_mapping(&p_lo, p_lo_row_mapper, q_lo_col_mapper, r_lo);
                inverse_mapping(&p_hi, p_hi_row_mapper, q_hi_col_mapper, r_hi);

                std::vector<int> good_row_pos;
                std::vector<int> good_col_pos;
                ant_passage(r_lo, r_hi, n, good_row_pos, good_col_pos);

                // merge r_hi to r_lo
                for (int i = 0; i < n; ++i) {
                    auto col = r_hi->get_col_by_row(i);
                    if (col == NOPOINT) continue;
                    r_lo->set_point(i, col);
                }

                // add good points
                for (int i = 0; i < good_col_pos.size(); ++i) {
                    auto col = good_col_pos[i];
                    auto row = good_row_pos[i];
                    r_lo->set_point(row, col);
                }


            }
        }


        /**
         * An implementation of staggered sticky multiplication aka glues two sticky braids to the new one  when
         * both have k common strands. Logic is as follows. We need to perform sticky multiplication on k common strands, while.
         * other parts should be untouched and remains the same. See details in the  book.
         * @param p
         * @param q
         * @param k
         * @param map
         * @param product
         */
        template<bool RowGlue, bool UseParallel>
        void staggered_sticky_multiplication(AbstractPermutation *p, AbstractPermutation *q, int k,
                                             PrecalcMap &map, AbstractPermutation *product,
                                             int nested_parall_regions = 0) {
            using namespace _details;


            if (k == p->row_size && k == q->col_size) {
                std::cout << "This function should not be called for this case, handled separately";
                return;
            }

            if (k == 0) {

                for (int i = 0; i < p->row_size; ++i) {
                    auto col = p->get_col_by_row(i);

                    if (RowGlue) {
                        if (col != NOPOINT) product->set_point(i + q->row_size, col + q->col_size);
                    } else {
                        if (col != NOPOINT) product->set_point(i, col);
                    }

                }

                for (int i = 0; i < q->row_size; ++i) {
                    auto col = q->get_col_by_row(i);

                    if (RowGlue) {
                        if (col != NOPOINT) product->set_point(i, col);
                    } else {
                        if (col != NOPOINT) product->set_point(i + p->row_size, col + p->col_size);
                    }

                }

                return;
            }


            int nearest_2_degree = pow(2, int(ceil(log2(2 * k))));
            int total = int(log2(nearest_2_degree)) * nearest_2_degree;

            int *memory_block;
            // then we need to use O(nlogn) memory
            if (UseParallel) {
                memory_block = new int[k * 8 + int(log2(nearest_2_degree)) * nearest_2_degree];
            } else {
                memory_block = new int[k * 8 + k * 8];
            }

            auto mapping_row = new int[k];
            auto mapping_col = new int[k];


            auto p_red = PermutationPreAllocated(k, k, memory_block, memory_block + k);
            auto q_red = PermutationPreAllocated(k, k, memory_block + 2 * k, memory_block + 3 * k);

            auto free_block_1 = memory_block;
            auto free_block_2 = memory_block + 4 * k;
            auto free_indices_block = memory_block + 8 * k;


            if (RowGlue) {
                // take first k columns from P and last k rows from Q, multiply and to bottom left corner of extended matrix
                get_vertical_slice(p, 0, k, mapping_row, &p_red);
                get_horizontal_slice(q, q->row_size - k, q->row_size, mapping_col, &q_red);
            } else {
                // take last k columns from P and first k rows from Q, multiply
                get_vertical_slice(p, p->row_size - k, p->row_size, mapping_row, &p_red);
                get_horizontal_slice(q, 0, k, mapping_col, &q_red);

            }

            if (UseParallel) {
                steady_ant_with_precalc_and_memory(&p_red, &q_red, free_block_1, free_block_2, free_indices_block, map,
                                                   total, nested_parall_regions);
            } else {
                steady_ant_optimized_seq(&p_red, &q_red, free_block_1, free_block_2, free_indices_block, map);
            }





            // result in p_red
            for (int i = 0; i < p_red.row_size; i++) {
                auto old_col = p_red.get_col_by_row(i);
                auto cur_col = mapping_col[old_col];
                auto cur_row = mapping_row[i];
                if (RowGlue) {
                    product->set_point(q->row_size - k + cur_row, cur_col);
                } else {
                    product->set_point(cur_row, cur_col + p->col_size - k);
                }
            }


            if (RowGlue) {
                for (int j = k; j < p->col_size; j++) {
                    auto row = p->get_row_by_col(j);
                    if (row != NOPOINT) product->set_point(row + q->row_size - k, j + q->col_size - k);
                }
            } else {
                for (int j = 0; j < p->col_size - k; j++) {
                    auto row = p->get_row_by_col(j);
                    if (row != NOPOINT) product->set_point(row, j);
                }
            }

            if (RowGlue) {
                for (int i = 0; i < q->row_size - k; i++) {
                    auto col = q->get_col_by_row(i);
                    if (col != NOPOINT) product->set_point(i, col);
                }
            } else {
                for (int i = k; i < q->row_size; i++) {
                    auto col = q->get_col_by_row(i);
                    if (col != NOPOINT) product->set_point(i - k + p->row_size, p->col_size - k + col);
                }
            }


            delete[] memory_block;
            delete[] mapping_col;
            delete[] mapping_row;


        }


        /**
         * An implementation of staggered sticky multiplication aka glues two sticky braids to the new one  when
         * both have k common strands. Logic is as follows. We need to perform sticky multiplication on k common strands, while.
         * other parts should be untouched and remains the same. See details in the  book.
         * @param p
         * @param q
         * @param k
         * @param map
        * @param product
        */
        void glueing_part_to_whole(AbstractPermutation *whole, AbstractPermutation *part, PrecalcMap &map, int offset_l,
                                   int offset_r, AbstractPermutation *product, int nested_parall_regions = 0) {
            using namespace _details;

            if (part->row_size != (whole->row_size - offset_l - offset_r)) {
                throw std::runtime_error("Dimensions not match");
            }

            auto k = whole->row_size - offset_r - offset_l;

            // first offset_l strands goes inact
            for (int col = 0; col < offset_l; ++col) {
                auto row = whole->get_row_by_col(col);
                product->set_point(row, col);
            }

            // last offset_r strands goes inact
            for (int col = whole->row_size - offset_r; col < whole->row_size; ++col) {
                auto row = whole->get_row_by_col(col);
                product->set_point(row, col);
            }


            int nearest_2_degree = pow(2, int(ceil(log2(2 * k))));
            int total = int(log2(nearest_2_degree)) * nearest_2_degree;
            auto memory_block = new int[k * 8 + int(log2(nearest_2_degree)) * nearest_2_degree];
            auto mapping_row = new int[k];

            auto whole_red = PermutationPreAllocated(k, k, memory_block, memory_block + k);
            auto part_copy = PermutationPreAllocated(k, k, memory_block + 2 * k, memory_block + 3 * k);
            whole_red.unset_all();
            part_copy.unset_all();

            copy(*part, part_copy);
            // get strands that insersects  with part braid  (its strands that hit border with part braid)
            get_vertical_slice(whole, offset_l, k + offset_l, mapping_row, &whole_red);


            steady_ant_with_precalc_and_memory(&whole_red, &part_copy, memory_block, memory_block + 4 * k,
                                               memory_block + 8 * k, map, total, nested_parall_regions);

            for (int i = 0; i < whole_red.row_size; i++) {
                auto old_col = whole_red.get_col_by_row(i);
                auto cur_col = old_col;
                auto cur_row = mapping_row[i];
                product->set_point(cur_row, cur_col + offset_l);
            }

            delete[] memory_block;
            delete[] mapping_row;
        }


        /**
         * Wrapper for steady ant multiplication, uses nlogn memory with since uses parallel approach
         * @param p permutation matrix that represents braid
         * @param q permutation matrix that represents braid
         * @param product product of two braids
         * @param map precalc values for all products of two braids up to specific size
         * @param nested_lvls number of upper levels of recursion  to parallelize. If 0 is set then algorithm follows sequential appoarch and uses linear space
         */
        void steady_ant_parallel_wrapper(AbstractPermutation &p, AbstractPermutation &q, AbstractPermutation &product,
                                         PrecalcMap &map, int nested_lvls = 3) {


            int nearest_2_degree = pow(2, int(ceil(log2(2 * p.row_size))));
            int total = int(log2(nearest_2_degree)) * nearest_2_degree;

            // if sequential mode
            if (nested_lvls == 0) total = p.row_size * 8;


            auto memory_block = new int[p.row_size * 8 + total];


            auto m_new = PermutationPreAllocated(p.row_size, p.col_size, memory_block, memory_block + p.row_size);
            auto n_new = PermutationPreAllocated(p.row_size, p.col_size, memory_block + 2 * p.row_size,
                                                 memory_block + 3 * p.row_size);

            copy(p, m_new);
            copy(q, n_new);

            if(nested_lvls == 0) {
                steady_ant_optimized_seq(
                        &m_new, &n_new, memory_block, memory_block + 4 * p.row_size, memory_block + 8 * p.row_size, map);
            } else {
                steady_ant_with_precalc_and_memory(
                        &m_new, &n_new, memory_block, memory_block + 4 * p.row_size, memory_block + 8 * p.row_size, map, total, nested_lvls);
            }

            copy(m_new, product);

            delete[] memory_block;

        }


    }


};

/**
  * For each size n (1,2,3...max_size) function compute
  * for all possible permutations p of size n  (n!) and all possible permutations of size n (n!)
  * it computes its distance product.
  * For example, map[3][hash_p][hash_q] will store product of permutation matrices p and q of size 3x3
  * @param map
  * @param max_size
  */
void precalc(PrecalcMap & map, int max_size) {
    using namespace std;
    int p_arr[max_size];
    int q_arr[max_size];
    auto empty_map = PrecalcMap();
    for (int size = 1; size < max_size + 1; size++) {
        for (int i = 0; i<size; ++i) p_arr[i] = i;


        do {
            for (int i = 0; i < size; ++i) q_arr[i] = i;
            do {
                auto p = new Permutation(size, size);
                for ( int i = 0; i < size; ++i) p->set_point(i, p_arr[i]);
                auto q = new Permutation(size, size);
                for (int i = 0;i < size; ++i) q-> set_point(i, q_arr[i]);

                long long hash_p = hash<AbstractPermutation>()(*p);
                long long hash_q = hash<AbstractPermutation>()(*q);
                auto r = distance_unit_monge_product::steady_ant::steady_ant(p, q, empty_map);
                auto points = std::vector<std::pair<int, int>>();
                if (map[size][hash_p].count(hash_q) > 0) {
                    std::cout << " Some error";
                    return;
                }
                r->to_points_on_grid(points);
                map[size][hash_p][hash_q] = points;
                delete r;
            }
            while (std::next_permutation(q_arr, q_arr + size));

        } while (
                std::next_permutation(p_arr, p_arr + size));
    }
}


#endif //CPU_STEADY_ANT_H
