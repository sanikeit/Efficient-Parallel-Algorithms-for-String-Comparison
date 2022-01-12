//
// File contains implementation of dominance sum incremental counting queries.
// Given permutation matrix with position with value at some position in associated monge matrix (dominance sum)
// functions calculate value in adjacent cells in O(1) time.
//

#include "matrices.h"

#ifndef CPU_DOMINANCE_SUMS_QUERIES_H
#define CPU_DOMINANCE_SUMS_QUERIES_H

#endif //CPU_DOMINANCE_SUMS_QUERIES_H

/**
 * Namespace contains functions that implements different incremental dominance sum queries
 */
namespace dominance_sum_counting {

    namespace top_left_arrow {

        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum++;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap >= col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap >= col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }


    }
    //ok
    namespace bottom_right_arrow {

        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);

            if (row_cap != NOPOINT && row_cap < row) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap < row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

    namespace top_right_arrow {


        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);

            if (row_cap >= row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

};
