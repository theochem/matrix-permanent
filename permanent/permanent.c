#include <stdbool.h>
#include <stdint.h>
#include <tgmath.h>

#include "permanent.h"
#include "binom.h"

#ifdef TUNING_FILE
/* Include tuning file. */
#include "tuning.h"
#else
/* Set default tuning parameters. */
#define PARAM_1 3.1415926
#define PARAM_2 8192
#endif


inline void swap2(int64_t *perm, int64_t i, int64_t j)
{
    const int64_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}


void init_perm(const int64_t N, int64_t *const the_fact_set, int64_t *const the_perm_set, int64_t *const the_inv_set)
{
    for (int64_t i = 0; i < N; i++) {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}


bool gen_next_perm(int64_t *const falling_fact, int64_t *const perm_, int64_t *const inv_perm, const int64_t cols, const int64_t u_)
{
    /* Use the current permutation to check for next permutation
     * lexicographically, if it exists update the curr_array by swapping the
     * leftmost changed element (at position i < k). Replace the elements up to
     * position k by the smallest element lying to the right of i. Do not put
     * remaining k, .., n - 1 positions in ascending order to improve efficiency
     * has time complexity O(n) in the worst case. */
    int64_t i = u_ - 1;
    int64_t m1 = cols - i - 1;

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
    int64_t z = (int64_t)perm_[i];
    do {
        ++z;
    } while (inv_perm[z] <= i);
    const int64_t j = inv_perm[z];
    swap2(perm_, i, j);
    swap2(inv_perm, perm_[i], perm_[j]);
    ++i;
    int64_t y = 0;

    /* Find the smallest elements to the right of position i. */
    while (i < u_) {
        while (inv_perm[y] < i) {
            ++y;
        }
        const int64_t j = inv_perm[y];
        swap2(perm_, i, j);
        swap2(inv_perm, perm_[i], perm_[j]);
        ++i;
    }
    return true;
}


double opt(const int64_t m_rows, const int64_t n_cols, const double *ptr)
{
    /* Use the fastest algorithm. */
    /* NOTE: This is just nonsense to show how it's done... */
    if (m_rows < PARAM_1) {
        return combinatoric(m_rows, n_cols, ptr);
    } else if (m_rows * n_cols < PARAM_2) {
        return glynn(m_rows, n_cols, ptr);
    } else {
        return ryser(m_rows, n_cols, ptr);
    }
}


double combinatoric(const int64_t m_rows, const int64_t n_cols, const double *ptr)
{
    /* sort up to position u + 1 where u = min(k, n_cols - 1). */
    int64_t sort_up_to = n_cols - 1;
    if (m_rows < n_cols) {
        sort_up_to = m_rows;
    }
    double sum_permanent = 0.0;
    double prod_permanent = 1.0;

    /* Allocate falling_fact, perm_, and inv_perm arrays. */
    int64_t falling_fact[128];
    int64_t perm_[128];
    int64_t inv_perm[128];
    /* Set sentinal value. */
    falling_fact[0] = 0;

    /* Initialize the set to permute. */
    init_perm(n_cols, falling_fact, perm_, inv_perm);
    bool gen_next_perm();

    /* Handle first permutation. */
    for (int64_t i = 0; i < m_rows; ++i) {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
    }
    sum_permanent = prod_permanent;

    /* Iterate over second to last permutations. */
    while (gen_next_perm(falling_fact, perm_, inv_perm, n_cols, sort_up_to)) {
        prod_permanent = 1.0;
        for (int64_t i = 0; i < m_rows; i++) {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
        sum_permanent += prod_permanent;
    }
    return sum_permanent;
}


double glynn(const int64_t m_rows, const int64_t n_cols, const double *ptr)
{
    /* Initialize gray code. */
    int64_t pos = 0;
    int64_t sign = 1;
    int64_t bound = n_cols - 1;
    int64_t delta[128];
    int64_t gray[128];

    /* Allocate and fill delta array (all +1 to start), and gray array from 0 to n_cols. */
    for (int64_t i = 0; i < n_cols; i++) {
        delta[i] = 1;
        gray[i] = i;
    }

    /* Allocate matrix for result of manual multiplies. */
    double vec[128];

    /* Dealing with a square matrix */
    if (m_rows == n_cols) {
        /* Handle first Gray code. */
        double result = 1.0;
        for (int64_t j = 0; j < n_cols; j++) {
            double sum_ = 0.0;
            for (int64_t i = 0; i < n_cols; i++) {
                /* Sum over all the values in each column. */
                sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
            }
            vec[j] = sum_;
        }
        for (int64_t j = 0; j < n_cols; j++) {
            result *= vec[j];
        }

        /* Iterate over second to last Gray codes. */
        while (pos != bound) {
            /* Update sign and delta. */
            sign *= -1;
            *(delta + bound - pos) *= -1;
            /* Compute each Gray code term. */
            for (int64_t j = 0; j < n_cols; j++) {
                double sum_ = 0.0;
                for (int64_t i = 0; i < n_cols; i++) {
                    sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
                }
                vec[j] = sum_;
            }
            double prod = 1.0;
            for (int64_t i = 0; i < n_cols; i++) {
                prod *= vec[i];
            }

            /* Multiply by the product of the vectors in delta. */
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

    /* Dealing with a rectangle. */
    else {
        /* Handle first Gray code. */
        double result = 1.0;
        for (int64_t j = 0; j < n_cols; j++) {
            double sum_ = 0.0;
            for (int64_t i = 0; i < m_rows; i++) {
                sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
            }
            for (int64_t k = m_rows; k < n_cols; k++) {
                sum_ += (double)delta[k];
            }
            vec[j] = sum_;
        }
        for (int64_t i = 0; i < n_cols; i++) {
            result *= vec[i];
        }

        /* Iterate over second to last Gray codes. */
        while (pos != bound) {
            /* Update sign and delta. */
            sign *= -1;
            *(delta + bound - pos) *= -1;
            /* Compute each Gray code term. */
            for (int64_t j = 0; j < n_cols; j++) {
                double sum_ = 0.0;
                for (int64_t i = 0; i < m_rows; i++) {
                    sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
                }

                for (int64_t k = m_rows; k < n_cols; k++) {
                    sum_ += (double)delta[k];
                }
                vec[j] = sum_;
            }
            double prod = 1.0;
            for (int64_t i = 0; i < n_cols; i++) {
                prod *= vec[i];
            }
            result += sign * prod;

            /* Go to next Gray code. */
            *gray = 0;
            *(gray + pos) = *(gray + pos + 1);
            ++pos;
            *(gray + pos) = pos;
            pos = gray[0];
        }

        /* Divide by external factor and return permanent. */
        return result / (pow(2.0, (double)bound) * (double)tgamma(n_cols - m_rows + 1));
    }
}


double ryser(const int64_t m_rows, const int64_t n_cols, const double *ptr)
{
    /* Initialize all relevant variables. See combinatorial algorithm for more
     * details as it was already went over. */
    int64_t falling_fact[128];
    int64_t perm_[128];
    int64_t inv_perm[128];
    double vec[128];
    falling_fact[0] = 0;

    /* Dealing with a square matrix. This bit-hacking trick was modified from
     * C++ code from Michelle Richer (lines 393-428) */
    if (m_rows == n_cols) {
        int32_t i, j, k;
        int64_t sum = 0, rowsum, rowsumprod;

        /* Iterate over c = pow(2, n) submatrices (equal to (1 << n)) submatrices. */
        int32_t c = 1UL << n_cols;

        /* Loop over columns of submatrix; compute product of row sums. */
        for (k = 0; k < c; k++) {
            rowsumprod = 1;
            for (i = 0; i < n_cols; i++) {
                /* Loop over rows of submatrix; compute row sum. */
                rowsum = 0;
                for (j = 0; j < n_cols; j++) {
                    /* Add element to row sum if the row index is in the characteristic vector of the submatrix, which is the binary vector given by k. */
                    if (k & (1UL << j)) {
                        rowsum += ptr[n_cols * i + j];
                    }
                }
                /* Update product of row sums. */
                rowsumprod *= rowsum;
            }
            /* Add term multiplied by the parity of the characteristic vector. */
            sum += rowsumprod * (1 - ((__builtin_popcount(k) & 1) << 1));
        }
        /* Return answer with the correct sign (times -1 for odd n). */
        int32_t sign = 1;
        if (n_cols % 2 == 1) {
            /* n is odd. */
            sign *= -1;
        }
        return sum * sign;
    }

    /* Dealing with a rectangle. Can't use bit hacking trick here. */
    else {
        int32_t value_sign = 1;
        int32_t counter = 0; // Count how many permutations you have generated
        double sum_over_k_vals = 0.0;
        for (int64_t k = 0; k < m_rows; k++) {
            /* Store the binomial coefficient for this k value bin_c. */
            double bin_c = BINOM(n_cols - m_rows + k, k);
            counter = 0;
            double sum_of_matrix_vals = 0.0;
            double prod_of_cols = 1.0;
            double result = 0.0;

            /* (Re)initialize the set to permute for this k value. */
            init_perm(n_cols, falling_fact, perm_, inv_perm);
            bool gen_next_perm();

            /* sort up to position u + 1 where u = min(m_rows - k, n_cols - 1). */
            int64_t sort_up_to = n_cols - 1;

            if ((m_rows - k) < sort_up_to) {
                sort_up_to = (m_rows - k);
            }

            for (int64_t i = 0; i < m_rows; i++) {
                sum_of_matrix_vals = 0.0;
                for (int64_t j = 0; j < sort_up_to; j++) {
                    sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
                }
                vec[i] = sum_of_matrix_vals;
                sum_of_matrix_vals = 0.0;
            }
            for (int64_t i = 0; i < m_rows; i++) {
                prod_of_cols *= vec[i];
            }

            result += value_sign * (double)bin_c * prod_of_cols;
            counter += 1;

            /* Iterate over second to last permutations of the set. */
            while (gen_next_perm(falling_fact, perm_, inv_perm, n_cols, sort_up_to)) {
                sum_of_matrix_vals = 0.0;
                for (int64_t i = 0; i < m_rows; i++) {
                    sum_of_matrix_vals = 0.0;
                    for (int64_t j = 0; j < m_rows - k; j++) {
                        sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
                    }
                    vec[i] = sum_of_matrix_vals;
                    sum_of_matrix_vals = 0.0;
                }
                prod_of_cols = 1.0;
                for (int64_t i = 0; i < m_rows; i++) {
                    prod_of_cols *= vec[i];
                }

                result += value_sign * (double)bin_c * prod_of_cols;
                counter += 1;
            }
            sum_over_k_vals += result / (counter / BINOM(n_cols, m_rows - k));
            value_sign *= -1;
        }

        return sum_over_k_vals;
    }
}
