#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "permanent.h"
#include "binomial.h"

#ifdef TUNING_FILE
/* Include tuning file. */
#include "tuning.h"
#else
/* Set default tuning parameters. */
#define PARAM_1 3.1415926
#define PARAM_2 8192
#endif


inline void swap2(size_t *, size_t, size_t);


inline void init_perm(const size_t, size_t *const, size_t *const, size_t *const);


bool gen_next_perm(size_t *const, size_t *const, size_t *const, const size_t, const size_t);


double opt(const int m_rows, const int n_cols, double *const ptr)
{
    /* Use the fastest algorithm. */
    if (m_rows < PARAM_1) {
        return combinatoric(m_rows, n_cols, ptr);
    } else if (m_rows * n_cols < PARAM_2) {
        return glynn(m_rows, n_cols, ptr);
    } else {
        return ryser(m_rows, n_cols, ptr);
    }
}


double combinatoric(const int m_rows, const int n_cols, double *const ptr)
{
    return (m_rows == n_cols) ? glynn(m_rows, n_cols, ptr) : glynn_rectangle(m_rows, n_cols, ptr);
}


double ryser(const int m_rows, const int n_cols, double *const ptr)
{
    return glynn(m_rows, n_cols, ptr);
}


double ryser_rectangle(const int m_rows, const int n_cols, double *const ptr)
{
    return glynn_rectangle(m_rows, n_cols, ptr);
}


double glynn(const int m_rows, const int n_cols, double *const ptr)
{
    size_t i, j;
    size_t pos = 0;
    size_t bound = n_cols - 1;
    size_t gray[128];

    int sign = 1;
    int delta[128];

    double sum, prod;
    double result = 1.0;
    double vec[128];

    /* Fill delta array (all +1 to start), and gray array from 0 to n_cols. */

    for (i = 0; i < n_cols; ++i) {
        delta[i] = 1;
        gray[i] = i;
    }

    /* Handle first Gray code. */

    for (j = 0; j < n_cols; ++j) {
        sum = 0.0;
        for (i = 0; i < n_cols; ++i) {
            /* Sum over all the values in each column. */
            sum += ptr[i * n_cols + j] * delta[i];
        }
        vec[j] = sum;
    }
    for (j = 0; j < n_cols; ++j) {
        result *= vec[j];
    }

    /* Iterate from the second to the final Gray code. */

    while (pos != bound) {
        /* Update sign and delta. */
        sign *= -1;
        delta[bound - pos] *= -1;
        // *(delta + bound - pos) *= -1;
        /* Compute each Gray code term. */
        for (j = 0; j < n_cols; ++j) {
            sum = 0.0;
            for (i = 0; i < n_cols; ++i) {
                sum += ptr[i * n_cols + j] * delta[i];
            }
            vec[j] = sum;
        }
        /* Multiply by the product of the vectors in delta. */
        prod = 1.0;
        for (i = 0; i < n_cols; ++i) {
            prod *= vec[i];
        }
        result += sign * prod;
        /* Go to next Gray code. */
        gray[0] = 0;
        gray[pos] = gray[pos + 1];
        ++pos;
        gray[pos] = pos;
        pos = gray[0];
    }

    /* Divide by external factor and return permanent. */

    return result / pow(2.0, (double)bound);
}


double glynn_rectangle(const int m_rows, const int n_cols, double *const ptr)
{
    size_t i, j, k;
    size_t pos = 0;
    size_t bound = n_cols - 1;
    size_t gray[128];

    int sign = 1;
    int delta[128];

    double sum, prod;
    double result = 1.0;
    double vec[128];

    /* Allocate and fill delta array (all +1 to start), and gray array from 0 to n_cols. */

    for (i = 0; i < n_cols; ++i) {
        delta[i] = 1;
        gray[i] = i;
    }

    /* Handle first Gray code. */

    for (j = 0; j < n_cols; ++j) {
        sum = 0.0;
        for (i = 0; i < m_rows; ++i) {
            sum += ptr[i * n_cols + j] * delta[i];
        }
        for (k = m_rows; k < n_cols; ++k) {
            sum += delta[k];
        }
        vec[j] = sum;
    }
    for (i = 0; i < n_cols; ++i) {
        result *= vec[i];
    }

    /* Iterate over second to last Gray codes. */

    while (pos != bound) {
        /* Update sign and delta. */
        sign *= -1;
        delta[bound - pos] *= -1;
        // *(delta + bound - pos) *= -1;
        /* Compute each Gray code term. */
        for (j = 0; j < n_cols; ++j) {
            sum = 0.0;
            for (i = 0; i < m_rows; ++i) {
                sum += ptr[i * n_cols + j] * delta[i];
            }

            for (k = m_rows; k < n_cols; ++k) {
                sum += delta[k];
            }
            vec[j] = sum;
        }
        /* Multiply by the product of the vectors in delta. */
        prod = 1.0;
        for (i = 0; i < n_cols; ++i) {
            prod *= vec[i];
        }
        result += sign * prod;
        /* Go to next Gray code. */
        gray[0] = 0;
        gray[pos] = gray[pos + 1];
        ++pos;
        gray[pos] = pos;
        pos = gray[0];
    }

    /* Divide by external factor and return permanent. */

    return result / (pow(2.0, (double)bound) * (double)tgamma(n_cols - m_rows + 1));
}


// double ryser(const int m_rows, const int n_cols, double *const ptr)
// {
//     /* Dealing with a square matrix. This bit-hacking trick was modified from
//      * C++ code from Michelle Richer (lines 289-316) */
//     size_t i, j;
//     unsigned long k;
//     double rowsum, rowsumprod;
//     double sum = 0.0;
//
//     /* Iterate over c = pow(2, n) submatrices (equal to (1 << n)) submatrices. */
//     unsigned long c = 1UL << n_cols;
//
//     /* Loop over columns of submatrix; compute product of row sums. */
//     for (k = 0; k < c; ++k) {
//         rowsumprod = 1.0;
//         for (i = 0; i < n_cols; ++i) {
//             /* Loop over rows of submatrix; compute row sum. */
//             rowsum = 0.0;
//             for (j = 0; j < n_cols; ++j) {
//                 /* Add element to row sum if the row index is in the characteristic vector of the submatrix, which is the binary vector given by k. */
//                 if (k & (1UL << j)) {
//                     rowsum += ptr[n_cols * i + j];
//                 }
//             }
//             /* Update product of row sums. */
//             rowsumprod *= rowsum;
//         }
//         /* Add term multiplied by the parity of the characteristic vector. */
//         sum += rowsumprod * (1 - ((__builtin_popcount(k) & 1) << 1));
//     }
//     /* Return answer with the correct sign (times -1 for odd n). */
//     return sum * ((n_cols % 2 == 1) ? 1 : -1);
// }
//
//
// double ryser_rectangle(const int m_rows, const int n_cols, double *const ptr)
// {
//     /* Initialize all relevant variables. See combinatorial algorithm for more
//      * details as it was already went over. */
//     size_t i, j, k, counter;
//     size_t sort_up_to;
//     size_t falling_fact[128];
//     size_t perm[128];
//     size_t inv_perm[128];
//     falling_fact[0] = 0;
//
//     int value_sign = 1;
//
//     double bin_c, sum_of_matrix_vals, prod_of_cols, result;
//     double sum_over_k_vals = 0.0;
//     double vec[128];
//
//     /* Dealing with a rectangle. Can't use bit hacking trick here. */
//     for (k = 0; k < m_rows; ++k) {
//         /* Store the binomial coefficient for this k value bin_c. */
//         bin_c = (double)(BINOM((n_cols - m_rows + k), k));
//         counter = 0; // Count how many permutations you have generated
//         sum_of_matrix_vals = 0.0;
//         prod_of_cols = 1.0;
//         result = 0.0;
//
//         /* (Re)initialize the set to permute for this k value. */
//         init_perm(n_cols, falling_fact, perm, inv_perm);
//
//         /* sort up to position u + 1 where u = min(m_rows - k, n_cols - 1). */
//         sort_up_to = n_cols - 1;
//
//         if ((m_rows - k) < sort_up_to) {
//             sort_up_to = (m_rows - k);
//         }
//
//         for (i = 0; i < m_rows; ++i) {
//             sum_of_matrix_vals = 0.0;
//             for (j = 0; j < sort_up_to; ++j) {
//                 sum_of_matrix_vals += ptr[i * n_cols + perm[j]];
//             }
//             vec[i] = sum_of_matrix_vals;
//             sum_of_matrix_vals = 0.0;
//         }
//         for (i = 0; i < m_rows; ++i) {
//             prod_of_cols *= vec[i];
//         }
//
//         result += value_sign * bin_c * prod_of_cols;
//         counter += 1;
//
//         /* Iterate over second to last permutations of the set. */
//         while (gen_next_perm(falling_fact, perm, inv_perm, n_cols, sort_up_to)) {
//             sum_of_matrix_vals = 0.0;
//             for (i = 0; i < m_rows; ++i) {
//                 sum_of_matrix_vals = 0.0;
//                 for (j = 0; j < m_rows - k; ++j) {
//                     sum_of_matrix_vals += ptr[i * n_cols + perm[j]];
//                 }
//                 vec[i] = sum_of_matrix_vals;
//                 sum_of_matrix_vals = 0.0;
//             }
//             prod_of_cols = 1.0;
//             for (i = 0; i < m_rows; ++i) {
//                 prod_of_cols *= vec[i];
//             }
//
//             result += value_sign * bin_c * prod_of_cols;
//             counter += 1;
//         }
//         sum_over_k_vals += (result / (double)counter) * BINOM(n_cols, (m_rows - k));
//         value_sign *= -1;
//     }
//     return sum_over_k_vals;
// }



inline void swap2(size_t *perm, size_t i, size_t j)
{
    const int temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}


inline void init_perm(const size_t N, size_t *const fact_set, size_t *const perm_set, size_t *const inv_set)
{
    size_t i;
    for (i = 0; i < N; i++) {
        fact_set[i + 1] = 0;
        perm_set[i] = i;
        inv_set[i] = i;
    }
}


bool gen_next_perm(size_t *const falling_fact, size_t *const perm, size_t *const inv_perm, const size_t cols, const size_t u)
{
    /* Use the current permutation to check for next permutation
     * lexicographically, if it exists update the curr_array by swapping the
     * leftmost changed element (at position i < k). Replace the elements up to
     * position k by the smallest element lying to the right of i. Do not put
     * remaining k, .., n - 1 positions in ascending order to improve efficiency
     * has time complexity O(n) in the worst case. */
    int i = u - 1;
    int m1 = cols - i - 1;

    /* Begin update of falling_fact - recall falling_fact[0] = 0, so increment
     * index accordingly. If i becomes negative during the check, you are
     * pointing to the sentinel value, so you are done generating
     * permutations. */
    while (falling_fact[i + 1] == m1) {
        falling_fact[i + 1] = 0;
        ++m1;
        --i;
    }
    if (i == -1) {
        return false;
    }
    ++falling_fact[i + 1];

    /* Find smallest element perm_[i] < perm_[j] that lies to the right of
     * pos i, and then update the state of the permutation using its inverse
     * to generate next. */
    int z = perm[i];
    do {
        ++z;
    } while (inv_perm[z] <= i);
    const int j = inv_perm[z];
    swap2(perm, i, j);
    swap2(inv_perm, perm[i], perm[j]);
    ++i;
    int y = 0;

    /* Find the smallest elements to the right of position i. */
    while (i < u) {
        while (inv_perm[y] < i) {
            ++y;
        }
        const int j = inv_perm[y];
        swap2(perm, i, j);
        swap2(inv_perm, perm[i], perm[j]);
        ++i;
    }
    return true;
}
