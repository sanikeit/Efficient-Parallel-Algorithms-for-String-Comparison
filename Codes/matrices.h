// File contains definition and implementation of permutation matrices
//
//

#ifndef CPU_MATRICES_H
#define CPU_MATRICES_H


#include <vector>

// NOPOINT indicates that no non-zero element in matrix in specified position
#define NOPOINT (-1)


/**
 * An abstract class that represents arbitrary matrix
 */
class AbstractMatrix {

public:
    int row_size;
    int col_size;

    virtual inline int get_element_at(int row, int col) const = 0;

    void print(std::ostream &os) const {

        for (int i = 0; i < row_size; ++i) {
            for (int j = 0; j < col_size; ++j) {
                os << get_element_at(i, j);
            }
            os << std::endl;
        }
    }

protected:
    AbstractMatrix(int row_size, int col_size) : row_size(row_size), col_size(col_size) {};

    virtual ~AbstractMatrix() = default;
};


/**
 * Class that represents permutation and sub-permutations matrices.
 * (sub-)permutation matrix is the matrix where in each row and each column there is (at most) exact one non-zero element
 * an example of permutation matrix: 0 1
 *                                   1 0
 * an example of subpermutation matrix:
 *                                   0 1
 *                                   0 0
 * NOPOINT indicates that no non-zero element in matrix in specified position
 */
class AbstractPermutation : public AbstractMatrix {
public:

    //TODO bad design but we need it for several algos
    int m = 0;
    int n = 0;

    /**
     * set non-zero element in (row,col) position in permutation matrix
     * @param row
     * @param col
     */
    inline virtual void set_point(int row, int col) = 0;


    /**
     * set value in position (row,col) to zero value
     * @param row
     * @param col
     */
    inline virtual void unset_point(int row, int col) = 0;


    /**
     * set values of all points of matrix to zero
     */
    inline virtual void unset_all() = 0;


    /**
     * if the matrix A contains a non-zero entry in the given column then the method returns
     * the associated row number with this non-zero entry such that A[row, col]=1, NOPOINT otherwise
     * @param col
     * @return row number of non-zero entry or NOPOINT
     */
    inline virtual int get_row_by_col(int col) const = 0;

    /**
    * if the matrix A contains a non-zero entry in the given row then the method returns
    * the associated col number with this non-zero entry such that A[row, col]=1 or NOPOINT otherwise
    * @param row
    * @return col number of non-zero entry or NOPOINT
     */
    inline virtual int get_col_by_row(int row) const = 0;

    virtual ~AbstractPermutation() = default;;


    AbstractPermutation(int row, int col) : AbstractMatrix(row, col) {}

    int get_element_at(int row, int col) const override {
        return (get_col_by_row(row) == col) ? 1 : 0;
    }


    /**
    * Checks whether or not this matrix equals to other matrix
    * @param other AbstractPermutation matrix
    * @return true if all non-zero entries in both matrices have same coordinates else false
    */
    bool is_equal_to(const AbstractPermutation &other) const {
        if (other.row_size != row_size || other.col_size != col_size) return false;
        for (int i = 0; i < row_size; ++i) {
            if (get_col_by_row(i) != other.get_col_by_row(i)) return false;
        }
        for (int i = 0; i < col_size; ++i) {
            if (get_row_by_col(i) != other.get_row_by_col(i)) return false;
        }
        return true;
    }


    /**
    * Fills input vector with position pairs of non-zero entries in the current matrix aka (row_i,col_i)
    * @param result std::vector<std::pair<int,int>>
    * @return void
    */
    void to_points_on_grid(std::vector<std::pair<int, int>> &result) const {
        for (int i = 0; i < col_size; ++i) {
            auto col = get_col_by_row(i);
            if (col != NOPOINT) result.emplace_back(i, col);
        }
    }

};


/**
 * Implementation of permutation matrix is based on two arrays:
 * the  first array is mapping of non-zero entries in rows to its position in cols
 * the second array is mapping of non-zero entries in cols to its position in rows
 * Note that memory management is handled inside class
 */
class Permutation : public AbstractPermutation {
private:
    int *row_to_col;
    int *col_to_row;

public:


    inline void set_point(int row, int col) override {
        row_to_col[row] = col;
        col_to_row[col] = row;
    }

    inline void unset_point(int row, int col) override {
        if (row_to_col[col] == col) {
            row_to_col[col] = NOPOINT;
            col_to_row[row] = NOPOINT;
        }
    }

    inline void unset_all() override {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }

    inline int get_row_by_col(int col) const override { return col_to_row[col]; }

    inline int get_col_by_row(int row) const override { return row_to_col[row]; }


    /**
     * Initializes permutation matrix of size rowXcol  as zero matrix that have no non-zero entries
     * @param row
     * @param col
     */
    explicit Permutation(int row, int col) : AbstractPermutation(row, col) {
        row_to_col = new int[row];
        col_to_row = new int[col];
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }


    /**
     * Initializes permutation matrix of size rowXcol  as matrix that has non-zero entries given by vector points
     * @param row
     * @param col
     * @param points
     */
    explicit Permutation(int row, int col, std::vector<std::pair<int, int>> &points) : Permutation(row, col) {
        for (auto &point: points) {
            row_to_col[point.first] = point.second;
            col_to_row[point.second] = point.first;
        }
    }

    ~Permutation() override {
        delete[] row_to_col;
        delete[] col_to_row;
    };


};


/**
 * Implementation of a permutation matrix is based on two arrays:
 * the  first array is mapping of non-zero entries in rows to its position in cols
 * the second array is mapping of non-zero entries in cols to its position in rows
 * Note that memory management is handled outside of  class. User provides memory block on which class is operates.
 * Memory freed outised of class.
 */
class PermutationPreAllocated : public AbstractPermutation {

public:

    int *row_to_col;
    int *col_to_row;

    inline void set_point(int row, int col) override {
        row_to_col[row] = col;
        col_to_row[col] = row;
    }

    inline void unset_point(int row, int col) override {
        if (row_to_col[col] == col) {
            row_to_col[col] = NOPOINT;
            col_to_row[row] = NOPOINT;
        }
    }

    inline void unset_all() override {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }


    inline int get_row_by_col(int col) const { return col_to_row[col]; }

    inline int get_col_by_row(int row) const { return row_to_col[row]; }


    PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row) : AbstractPermutation(row, col),
                                                                                  row_to_col(row_to_col),
                                                                                  col_to_row(col_to_row) {};


    PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row,
                            std::vector<std::pair<int, int>> &points) :
            PermutationPreAllocated(row, col, row_to_col, col_to_row) {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
        for (auto &point: points) {
            row_to_col[point.first] = point.second;
            col_to_row[point.second] = point.first;
        }
    }

    ~PermutationPreAllocated() override = default;;

};


/**
 * The class presents a simple 2d matrix with integer value entities. Used primarily for testing
 */
class Matrix : public AbstractMatrix {
private:
    int *arr;
public:

    Matrix(int row_size, int col_size) : AbstractMatrix(row_size, col_size) {
        arr = new int[row_size * col_size];
        for (int i = 0; i < row_size * col_size; ++i) arr[i] = 0;
    }

    int get_element_at(int row, int col) const override { return arr[row * col_size + col]; }

    inline void set_element_at(int row, int col, int value) { arr[row * col_size + col] = value; }


    void print(std::ostream &os) {

        for (int i = 0; i < row_size; ++i) {
            for (int j = 0; j < col_size; ++j) {
                os << get_element_at(i, j);
            }
            os << std::endl;
        }
    }

    ~Matrix() override {
        delete[] arr;
    }
};


/**
 * Fill given zero permutation matrix randomly by given seed
 * @param m
 * @param row_size
 * @param col_size
 * @param seed
 */
void fill_permutation_matrix(AbstractPermutation *m, int row_size, int col_size, int seed = 0) {

    auto is_used = new bool[col_size];
    for (int i = 0; i < col_size; ++i) is_used[i] = false;

    /* initialize random seed: */
    srand(seed);

    auto active_pts = col_size;

    for (int row = 0; row < row_size && active_pts > 0; ++row) {

        while (true) {
            auto col = abs(rand()) % col_size;
            if (!is_used[col]) {
                m->set_point(row, col);
                is_used[col] = true;
                break;
            }
        }
        active_pts--;
    }

    delete[] is_used;
}


/**
 * Get dominance matrix of specified func operator( could left/right bottom/top arrow) from permutations matrix m
 * @tparam Lambda
 * @param m
 * @param func
 * @param output_dominance_matrix
 */
template<typename Lambda>
void get_dominance_matrix(AbstractPermutation &m, Lambda &&func, Matrix *output_dominance_matrix) {
    auto row_size = m.row_size + 1;
    auto col_size = m.col_size + 1;

    for (int row = 0; row < row_size; ++row) {
        for (int col = 0; col < col_size; ++col) {
            for (int i = 0; i < m.row_size; ++i) {
                auto row_pos_point = i;
                auto col_pos_point = m.get_col_by_row(row_pos_point);
                if (col_pos_point == NOPOINT) continue;
                if (func(row_pos_point, col_pos_point, row, col) == true)
                    output_dominance_matrix->set_element_at(row, col,
                                                            output_dominance_matrix->get_element_at(row, col) + 1);
            }
        }
    }
}


/**
 * Cope elements from permutation matrix from to permutation matrix to
 * @param from
 * @param to
 */
inline void copy(AbstractPermutation &from, AbstractPermutation &to) {
    for (int i = 0; i < from.row_size; ++i) {
        auto col = from.get_col_by_row(i);
        to.set_point(i, col);
    }
}



#endif //CPU_MATRICES_H
