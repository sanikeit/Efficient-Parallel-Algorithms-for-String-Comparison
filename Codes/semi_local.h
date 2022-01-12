#ifndef CPU_SEMI_LOCAL_H
#define CPU_SEMI_LOCAL_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include <map>
#include <unordered_map>
#include "steady_ant.h"
#include  <cstdlib>
#include <chrono>
#include "predefined_types.h"
#include "deque"


namespace semi_local {

    namespace _details {

        /**
         *
         * @tparam Input
         * @tparam WithIf weather or not to use approach with if rather then  the branchless one @param strand_map
         * @param a
         * @param b
         * @param upper_bound
         * @param left_edge
         * @param top_edge
         * @param offset_a
         * @param offset_b
         */
        template<class Input, bool WithIf, bool withWait>
        inline void
        anti_diagonal_computation(Input *strand_map, const Input *a, const Input *b, int upper_bound, int left_edge,
                                  int top_edge, int offset_a, int offset_b) {

#pragma omp for simd schedule(static) aligned(a, b, strand_map:sizeof(Input)*8) nowait
            for (int k = 0; k < upper_bound; ++k) {

                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];

                auto r = (a[offset_a + k] == b[offset_b + k]) || (left_strand > right_strand);

                if (WithIf) {
                    if (r) {
                        strand_map[top_edge + k] = left_strand;
                        strand_map[left_edge + k] = right_strand;
                    }
                } else {
                    auto r_minus = (r - 1);
                    auto minus_r = -r;
                    auto l_new = (left_strand & r_minus) | (minus_r  & right_strand);
                    auto r_new = (right_strand & r_minus) | (minus_r & left_strand);

                    strand_map[left_edge + k] = l_new;
                    strand_map[top_edge + k] = r_new;
                }
            }

            if (withWait) {
#pragma omp barrier
            }
        }


        template<class Input>
        inline void initialization(Input *strand_map, int m, int n) {
#pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                strand_map[k] = k;
            }

#pragma omp for simd schedule(static)
            for (int l = 0; l < n; ++l) {
                strand_map[l + m] = l + m;
            }

        }

        template<class Input>
        inline void
        construct_permutation(AbstractPermutation &matrix, Input *strand_map, bool is_reverse, int m, int n) {
            if (!is_reverse) {
#pragma omp for simd schedule(static)
                for (int r = m; r < m + n; r++) {
                    matrix.set_point(strand_map[r], r - m);
                }
#pragma omp for simd schedule(static)
                for (int l = 0; l < m; l++) {
                    matrix.set_point(strand_map[l], n + l);
                }


            } else {

#pragma omp for simd schedule(static)
                for (int r = m; r < m + n; r++) {
                    matrix.set_point(n + m - 1 - strand_map[r], n + m - 1 - (r - m));
                }
#pragma omp for simd schedule(static)
                for (int l = 0; l < m; l++) {
                    matrix.set_point(n + m - 1 - strand_map[l], n + m - 1 - (n + l));
                }

            }

        }

        template<class Input>
        inline void fill_a_reverse(const Input *a, Input *a_reverse, int m) {
#pragma omp  for simd schedule(static)
            for (int i = 0; i < m; ++i) {
                a_reverse[i] = a[m - 1 - i];
            }
        }


    }


    /**
     * Computes the kernel of semi-local lcs solution for given two strings in naive fashion
     * @tparam Input
     * @tparam WithIf
     * @param permutation
     * @param a
     * @param a_size
     * @param b
     * @param b_size
     */
    template<class Input, bool WithIf>
    void
    sticky_braid_sequential(AbstractPermutation &permutation, const Input *a, int a_size, const Input *b, int b_size) {

        auto m = a_size;
        auto n = b_size;

        auto top_strands = new Input[n];
        auto left_strands = new Input[m];

        // init phase
        for (int i = 0; i < m; ++i) {
            left_strands[i] = i;
        }
        for (int i = 0; i < n; ++i) {
            top_strands[i] = i + m;
        }

        for (int i = 0; i < m; ++i) {
            auto left_edge = m - 1 - i;
            auto left_strand = left_strands[left_edge];
            auto a_symbol = a[i];
            int right_strand;
            for (int j = 0; j < n - 1; ++j) {
                right_strand = top_strands[j];
                auto r = a_symbol == b[j] || (left_strand > right_strand);

                if (WithIf) {
                    if (r) {
                        top_strands[j] = left_strand;
                        left_strand = right_strand;
                    }
                } else {

                    top_strands[j] = (right_strand & (r - 1)) | ((-r) & left_strand);
                    left_strand = (left_strand & (r - 1)) | ((-r) & right_strand);

                }

            }

            right_strand = top_strands[n - 1];
            auto r = a_symbol == b[n - 1] || (left_strand > right_strand);
            left_strands[left_edge] = (left_strand & (r - 1)) | ((-r) & right_strand);

            if (WithIf) {
                if (r) top_strands[n - 1] = left_strand;
            } else {
                top_strands[n - 1] = (right_strand & (r - 1)) | ((-r) & left_strand);
            }

        }


        // permutation construction phase
        for (int l = 0; l < m; l++) permutation.set_point(left_strands[l], n + l);

        for (int r = m; r < m + n; r++) permutation.set_point(top_strands[r - m], r - m);

        delete[] left_strands;
        delete[] top_strands;
    }


    /**
     * Computes the kernel of semi-local lcs solution for given two strings with antidiagonal pattern using Open MP
     * @tparam Input
     * @tparam WithIf
     * @tparam WithWait
     * @param matrix
     * @param a
     * @param a_size
     * @param b
     * @param b_size
     * @param threads_num
     * @param is_reverse
     */
    template<class Input, bool WithIf, bool WithWait>
    void sticky_braid_mpi(AbstractPermutation &matrix, const Input *a, int a_size, const Input *b, int b_size,
                          int threads_num = 1, bool is_reverse = false) {
        using namespace _details;
        using namespace distance_unit_monge_product::steady_ant::_details;

        if (a_size > b_size) {
            sticky_braid_mpi<Input, WithIf, WithWait>(matrix, b, b_size, a, a_size, threads_num, !is_reverse);
            return;
        }

        auto m = a_size;
        auto n = b_size;


        auto size = m + n;
        Input *strand_map = new Input[size];

        auto num_diag = m + n - 1;
        auto total_same_length_diag = num_diag - (m - 1) - (m - 1);
        Input *a_reverse = new Input[m];

#pragma omp parallel num_threads(threads_num)  default(none) shared(a_reverse, a, b, is_reverse, strand_map, matrix, total_same_length_diag, size, m, n)
        {


            int left_edge, top_edge;
            //    init phase
            initialization(strand_map, m, n);
            fill_a_reverse(a, a_reverse, m);

            //    phase one
            top_edge = m;
            left_edge = m - 1;
            for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
                anti_diagonal_computation<Input, WithIf, WithWait>(strand_map, a_reverse, b, cur_diag_len + 1,
                                                                   left_edge,
                                                                   top_edge, left_edge, 0);
                left_edge--;
            }

            //phase 2
            top_edge = m;
            for (int j = 0; j < total_same_length_diag; ++j) {
                anti_diagonal_computation<Input, WithIf, WithWait>(strand_map, a_reverse, b, m, 0, top_edge, 0, j);
                top_edge++;
            }

            //// phase 3
            auto start_j = total_same_length_diag;
            top_edge = start_j + m;

            for (int diag_len = m - 2; diag_len >= 0; --diag_len, start_j++) {
                anti_diagonal_computation<Input, WithIf, WithWait>(strand_map, a_reverse, b, diag_len + 1, 0, top_edge,
                                                                   0,
                                                                   start_j);
                top_edge++;
            }

            construct_permutation(matrix, strand_map, is_reverse, m, n);
        }

        delete[] a_reverse;
        delete[] strand_map;
    }

    /**
     * Computes the kernel of semi-local lcs solution for given two strings with antidiagonal pattern using Open MP where
     * 1st and 3rd phases are merged together so there less syncronizations
     * @tparam Input
     * @tparam WithIf
     * @param matrix
     * @param a
     * @param a_size
     * @param b
     * @param b_size
     * @param map
     * @param nested_parall_regions
     * @param threads_num
     */
    template<class Input, bool WithIf, bool WithWait>
    void
    first_and_third_phase_combined(AbstractPermutation &matrix, const Input *a, int a_size, const Input *b,
                                   int b_size, PrecalcMap &map, int nested_parall_regions = 0, int threads_num = 1) {
        using namespace distance_unit_monge_product::steady_ant::_details;
        using namespace _details;
        using distance_unit_monge_product::steady_ant::glueing_part_to_whole;

        if (a_size > b_size) {
            auto p = Permutation(a_size + b_size, a_size + b_size);
            first_and_third_phase_combined<Input, WithIf, WithWait>(p, b, b_size, a, a_size, map, nested_parall_regions,
                                                                    threads_num);
            fill_permutation_ba(&p, &matrix, a_size, b_size);
            return;
        }

        //assume |a|<=|b|

        auto m = a_size;
        auto n = b_size;

        auto size = m + n;
        Input *strand_map = new Input[size + 2 * (m - 1)];
        auto third_phase_map_size = m * 2 - 2;
        auto third_phase_map = strand_map + size;

        auto p = Permutation(m + n, m + n);
        auto q = Permutation(third_phase_map_size, third_phase_map_size);

        auto offset = n - (m - 1);
        Input *a_reverse = new Input[m];

#pragma omp parallel num_threads(threads_num)  default(none) shared(a, a_reverse, b, strand_map, size, m, n, matrix, p, q, offset, third_phase_map, third_phase_map_size)
        {

            fill_a_reverse(a, a_reverse, m);
            int in_third_phase = m - 1;
            //    init phase
#pragma omp for simd schedule(static) nowait
            for (int k = 0; k < (m + n); ++k) {
                strand_map[k] = k;
            }

#pragma omp for simd schedule(static) nowait
            for (int k = 0; k < third_phase_map_size; ++k) {
                if (k < m - 1) {
                    third_phase_map[k] = 2 * k;
                } else {
                    third_phase_map[k] = (k - (m - 1)) * 2 + 1;
                }
            }
#pragma omp barrier


            for (int diag_number = 0; diag_number < m - 1; ++diag_number) {


#pragma omp for simd schedule(static) nowait
                for (int pos_in_diag = 0; pos_in_diag < in_third_phase; ++pos_in_diag) {

                    auto top_edge = diag_number + pos_in_diag;
                    auto left_strand = third_phase_map[pos_in_diag];
                    auto top_strand = third_phase_map[m - 1 + top_edge];
                    bool r = a_reverse[pos_in_diag] == b[offset + top_edge] || (left_strand > top_strand);

                    if (WithIf) {
                        if (r) std::swap(third_phase_map[pos_in_diag], third_phase_map[m - 1 + top_edge]);
                    } else {
                        third_phase_map[pos_in_diag] = (left_strand & (r - 1)) | ((-r) & top_strand);
                        third_phase_map[m - 1 + top_edge] = (top_strand & (r - 1)) | ((-r) & left_strand);
                    }


                }

#pragma omp for simd schedule(static)
                for (int pos_in_diag = in_third_phase; pos_in_diag < m; ++pos_in_diag) {
                    auto top_edge = diag_number + pos_in_diag + 1 - m;
                    auto left_strand = strand_map[pos_in_diag];
                    auto top_strand = strand_map[m + top_edge];
                    bool r = a_reverse[pos_in_diag] == b[top_edge] || (left_strand > top_strand);

                    if (WithIf) {
                        if (r) if (r) std::swap(strand_map[pos_in_diag], strand_map[m + top_edge]);
                    } else {
                        strand_map[pos_in_diag] = (left_strand & (r - 1)) | ((-r) & top_strand);
                        strand_map[m + top_edge] = (top_strand & (r - 1)) | ((-r) & left_strand);
                    }


                }
                in_third_phase--;
            }

            //phase 2
            auto top_edge = m;
            for (int j = 0; j < offset; ++j) {
                anti_diagonal_computation<Input, WithIf, WithWait>(strand_map, a_reverse, b, m, 0, top_edge, 0, j);
                top_edge++;
            }

#pragma omp for simd schedule(static) nowait
            for (int l = 0; l < m; l++) {
                if (l == m - 1) {
                    p.set_point(strand_map[l], n + l);
                } else {
                    p.set_point(strand_map[l], l * 2 + offset);
                }
            }


#pragma omp for simd schedule(static) nowait
            for (int r = m; r < m + n; r++) {
                if ((r - m) < offset) {
                    p.set_point(strand_map[r], r - m);
                } else {
                    p.set_point(strand_map[r], (r - m - offset + 1) * 2 + offset - 1);
                }
            }

#pragma omp for simd schedule(static) nowait
            for (int l = 0; l < m - 1; l++) {
                q.set_point(third_phase_map[l], m - 1 + l);
            }


#pragma omp for simd schedule(static) nowait
            for (int r = m - 1; r < m + m - 2; r++) {
                q.set_point(third_phase_map[r], r - (m - 1));
            }

#pragma omp barrier
        }

        glueing_part_to_whole(&p, &q, map, offset, 1, &matrix, nested_parall_regions);

        delete[] strand_map;
        delete[] a_reverse;
    }


    /**
     * Hybrid appoarch of recursive and iterative combing. Firstly, follows the recursive structure, then switches to the iterative combing
     * TODO add documentation in details
     * @tparam Input
     * @tparam WithIf
     * @tparam WithWait
     * @tparam UseSumBound
     * @param perm
     * @param a
     * @param m
     * @param b
     * @param n
     * @param map
     * @param thds_per_combing_algo
     * @param braid_mul_parall_depth
     * @param depth
     * @param sum_bound
     * @param parallel_depth
     */
    template<class Input, bool WithIf, bool WithWait, bool UseSumBound>
    void hybrid(
            AbstractPermutation &perm, const Input *a, int m, const Input *b, int n,
            PrecalcMap &map, int thds_per_combing_algo,
            int braid_mul_parall_depth, int depth, int sum_bound, int parallel_depth) {

        using namespace _details;
        using distance_unit_monge_product::steady_ant::staggered_sticky_multiplication;


        if (UseSumBound) {
            if (m + n <= sum_bound) {
                sticky_braid_mpi<Input, WithIf, WithWait>(perm, a, m, b, n, thds_per_combing_algo);
                return;
            }

        } else {
            //base case
            if (depth <= 0) {
                sticky_braid_mpi<Input, WithIf, WithWait>(perm, a, m, b, n, thds_per_combing_algo);
                return;
            }
        }


        if (n > m) {
            auto n1 = n / 2;
            auto b1 = b;
            auto b2 = b + n1;

            auto subtree_l = Permutation(n1 + m, n1 + m);
            auto subtree_r = Permutation(n - n1 + m, n - n1 + m);

            if (parallel_depth > 0) {
#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {
#pragma omp task
                        hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_l, a, m, b1, n1, map,
                                                                     thds_per_combing_algo,
                                                                     braid_mul_parall_depth, depth - 1, sum_bound,
                                                                     parallel_depth - 1);

#pragma omp task
                        hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_r, a, m, b2, n - n1, map,
                                                                     thds_per_combing_algo,
                                                                     braid_mul_parall_depth, depth - 1, sum_bound,
                                                                     parallel_depth - 1);
                    }

                }
#pragma omp taskwait
            } else {
                hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_l, a, m, b1, n1, map, thds_per_combing_algo,
                                                             braid_mul_parall_depth,
                                                             depth - 1, sum_bound, parallel_depth - 1);
                hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_r, a, m, b2, n - n1, map, thds_per_combing_algo,
                                                             braid_mul_parall_depth, depth - 1, sum_bound,
                                                             parallel_depth - 1);
            }


            staggered_sticky_multiplication<false, false>(&subtree_l, &subtree_r, m, map, &perm,
                                                          braid_mul_parall_depth);

        } else {

            auto m1 = m / 2;
            auto a1 = a;
            auto a2 = a + m1;

            auto subtree_l = Permutation(m1 + n, m1 + n);
            auto subtree_r = Permutation(m - m1 + n, m - m1 + n);


            if (parallel_depth > 0) {
#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {
#pragma omp task
                        hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_l, a1, m1, b, n, map,
                                                                     thds_per_combing_algo,
                                                                     braid_mul_parall_depth, depth - 1, sum_bound,
                                                                     parallel_depth - 1);
#pragma omp task
                        hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_r, a2, m - m1, b, n, map,
                                                                     thds_per_combing_algo,
                                                                     braid_mul_parall_depth, depth - 1, sum_bound,
                                                                     parallel_depth - 1);
                    }
                }
#pragma omp taskwait
            } else {
                hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_l, a1, m1, b, n, map, thds_per_combing_algo,
                                                             braid_mul_parall_depth,
                                                             depth - 1, sum_bound, parallel_depth - 1);
                hybrid<Input, WithIf, WithWait, UseSumBound>(subtree_r, a2, m - m1, b, n, map, thds_per_combing_algo,
                                                             braid_mul_parall_depth, depth - 1, sum_bound,
                                                             parallel_depth - 1);
            }

            staggered_sticky_multiplication<true, false>(&subtree_l, &subtree_r, n, map, &perm, braid_mul_parall_depth);
        }


    }


    /**
     * The hybrid appaorch with down to top apparoch. No recustion. Fixed amount of iteartive combing problems that further are merged via
     * sticky braid multiplication.
     * @tparam Input
     * @tparam WithIf
     * @param perm
     * @param a
     * @param m
     * @param b
     * @param n
     * @param map
     * @param small_m
     * @param small_n
     * @param threads_num
     */
    template<class Input, bool WithIf>
    void hybrid_iterative(
            AbstractPermutation &perm, const Input *a, int m, const Input *b, int n,
            PrecalcMap &map, int small_m, int small_n, int threads_num = 1) {

        using namespace _details;
        using distance_unit_monge_product::steady_ant::staggered_sticky_multiplication;


        int cols_per_block = n / small_n;
        int rows_per_block = m / small_m;


        int num_tasks = small_m * small_n;

        auto tasks = new AbstractPermutation *[num_tasks];
        auto tasks_next_iter = new AbstractPermutation *[num_tasks];


#pragma omp parallel  master taskloop num_threads(threads_num)
        for (int i = 0; i < num_tasks; i++) {
            int start_col = (i % small_n) * cols_per_block;
            int end_col = std::min(start_col + cols_per_block, n);

            int start_row = (i / small_n) * rows_per_block;
            int end_row = std::min(start_row + rows_per_block, m);

            // the edge blocks may need do extra work
            if ((i % small_n) == small_n - 1) end_col = n;


            if ((i / small_n) == small_m - 1) end_row = m;


            int size_block_b = end_col - start_col;
            int size_block_a = end_row - start_row;


            auto b_loop = b + start_col;
            auto a_loop = a + start_row;

            auto matrix = new Permutation(size_block_a + size_block_b, size_block_a + size_block_b);

            matrix->m = size_block_a;
            matrix->n = size_block_b;

            sticky_braid_mpi<Input, WithIf, false>(*matrix, a_loop, size_block_a, b_loop, size_block_b, 1);


            tasks[i] = matrix;
        }


        auto next_jobs = tasks;
        auto current_jobs = tasks_next_iter;

        int steps = ceil(log2(small_n)) + ceil(log2(small_m));


        int block_rows = rows_per_block;
        int block_cols = cols_per_block;


        double total = 0;
        for (int j = 0; j < steps; ++j) {


            bool is_reduction_in_row = small_n > small_m;

            //TODO: heuristic now we choose to glue by the large one (seems to work) best. Specifically reduction in a row if on next step row is still big
            //TODO: following is working worse but it copipies logic of recursion: is_reduction_in_row = 2 * block_cols >= block_rows
            if (small_n > 1 && small_m > 1) {
                is_reduction_in_row = block_rows >= 2 * block_cols;
            }


            if (is_reduction_in_row) {
                block_cols *= 2;
            } else {
                block_rows *= 2;
            }


            auto new_cols = is_reduction_in_row ? int(ceil(small_n / 2.0)) : small_n;
            auto new_rows = !is_reduction_in_row ? int(ceil(small_m / 2.0)) : small_m;

            auto tmp = current_jobs;
            current_jobs = next_jobs;

            next_jobs = tmp;


#pragma omp parallel master taskloop num_threads(threads_num)
            for (int i = 0; i < new_cols * new_rows; i++) {

                auto cur_row = i / new_cols;
                auto cur_col = i % new_cols;


                AbstractPermutation *p;
                AbstractPermutation *q;

                if (is_reduction_in_row) {
                    p = current_jobs[cur_row * small_n + 2 * cur_col];

                    if ((2 * cur_col + 1) >= small_n) {
                        next_jobs[i] = p;
                    } else {
                        q = current_jobs[cur_row * small_n + 2 * cur_col + 1];
                        auto product = new Permutation(p->m + q->n + p->n, p->m + q->n + p->n);
                        product->m = p->m;
                        product->n = p->n + q->n;

                        staggered_sticky_multiplication<false, false>(p, q, p->m, map, product);

                        next_jobs[i] = product;

                        delete p;
                        delete q;
                    }

                } else {
                    p = current_jobs[2 * cur_row * small_n + cur_col];

                    if ((2 * cur_row + 1) >= small_m) {
                        next_jobs[i] = p;
                    } else {
                        q = current_jobs[(2 * cur_row + 1) * small_n + cur_col];
                        auto product = new Permutation(p->m + q->m + p->n, p->m + q->m + p->n);
                        product->m = p->m + q->m;
                        product->n = p->n;

                        staggered_sticky_multiplication<true, false>(p, q, p->n, map, product);

                        next_jobs[i] = product;

                        delete p;
                        delete q;
                    }


                }

            }

            small_n = new_cols;
            small_m = new_rows;

        }

        auto result = next_jobs[0];


        //        fill permutation perm
#pragma omp parallel for
        for (int j = 0; j < m + n; ++j) perm.set_point(j, result->get_col_by_row(j));

        delete result;
        delete[] tasks;
        delete[] tasks_next_iter;

    }

    template<bool WithIf>
    void
    hybrid_iterative_wrapper(AbstractPermutation &perm, const unsigned short *a, int m, const unsigned short *b, int n,
                             PrecalcMap &map, int threads_num = 1) {



        /**
         * Heuristic for spliting is as follows.
         * Set as small as possible blocks in one dimenstion bounded by either 32000 or m.
         * If m is greater then  small_m is  m / 32k
         *
         * For other dimenstion. We split equally work among threads. So they get same task with same length
         */

        int cols_per_block = 32000;
        int rows_per_block = 32000;

        int m_small;
        int n_small;

        // if m is less then bound then n_small would be 1
        if (rows_per_block >= m) {
            rows_per_block = m;
            m_small = 1;
        } else {
            // else we check how many blocks in row we need with size rows_per_block
            // and adjust rows_per_block to equal partion
            int nums = int(ceil((1.0 * m) / rows_per_block));
            rows_per_block = int(ceil(1.0 * m / nums));
            m_small = ceil((m * 1.0) / rows_per_block);
        }


        // if fits then equally split between threads
        if (cols_per_block >= n) {
            cols_per_block = int(ceil((n * 1.0) / threads_num));
            n_small = threads_num;
        } else {
            int cols_per_thread = int(ceil((1.0 * n) / threads_num));
            int rest = (64000 - rows_per_block);

            if (cols_per_thread < rest) {
                n_small = threads_num;
            } else {
                n_small = threads_num * ceil(1.0 * cols_per_thread / rest);
            }
        }

        hybrid_iterative<unsigned short, WithIf>(perm, a, m, b, n, map, m_small, n_small, threads_num);
    }


}

#endif //CPU_SEMI_LOCAL_H
