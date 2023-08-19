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


double opt(const size_t m_rows, const size_t n_cols, double *const ptr)
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

inline void swap2(size_t *perm, size_t i, size_t j)
{
    const size_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}


void init_perm(const size_t N, size_t *const the_fact_set, size_t *const the_perm_set, size_t *const the_inv_set)
{
    for (size_t i = 0; i < N; i++) {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}


bool gen_next_perm(size_t *const falling_fact, size_t *const perm, size_t *const inv_perm, const size_t cols, const size_t u_)
{
    /* Use the current permutation to check for next permutation
     * lexicographically, if it exists update the curr_array by swapping the
     * leftmost changed element (at position i < k). Replace the elements up to
     * position k by the smallest element lying to the right of i. Do not put
     * remaining k, .., n - 1 positions in ascending order to improve efficiency
     * has time complexity O(n) in the worst case. */
    size_t i = u_ - 1;
    size_t m1 = cols - i - 1;

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

    /* Find smallest element perm[i] < perm[j] that lies to the right of
     * pos i, and then update the state of the permutation using its inverse
     * to generate next. */
    size_t z = (size_t)perm[i];
    do {
        ++z;
    } while (inv_perm[z] <= i);
    const size_t j = inv_perm[z];
    swap2(perm, i, j);
    swap2(inv_perm, perm[i], perm[j]);
    ++i;
    size_t y = 0;

    /* Find the smallest elements to the right of position i. */
    while (i < u_) {
        while (inv_perm[y] < i) {
            ++y;
        }
        const size_t j = inv_perm[y];
        swap2(perm, i, j);
        swap2(inv_perm, perm[i], perm[j]);
        ++i;
    }
    return true;
}


double combinatoric(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    size_t i, j, tmp, bound;
    size_t gray[64 + 1];
    size_t perm[64];

    double prod;
    double result = 1.0;

    /* Fill permutation array and Gray code array with [0...n_cols]. */

    for (i = 0; i < n_cols; ++i) {
      gray[i] = i;
      perm[i] = i;
    }
    gray[n_cols] = n_cols;

    /* Handle first Gray code; compute first product of A_{i,\sigma(i)}. */

    for (i = 0; i < m_rows; ++i) {
      result *= ptr[i * n_cols + perm[i]];
    }

    /* Iterate from the second to the final Gray code. */

    i = 1;

    if (m_rows < n_cols)
        bound = m_rows + 1;
    else
        bound = n_cols;

    while (i < bound) {

        /* Advance to next permutation. */

        --gray[i];
        j = i % 2 * gray[i];
        tmp = perm[j];
        perm[j] = perm[i];
        perm[i] = tmp;

        /* Compute product of A_{i,\sigma(i)}. */

        prod = 1.0;
        for (j = 0; j < m_rows; j++) {
            prod *= ptr[j * n_cols + perm[j]];
        }
        result += prod;

        /* Go to next Gray code. */

        i = 1;
        while (!gray[i]) {
            gray[i] = i;
            ++i;
        }
    }

    return result;
}

double combinatoric_rectangle(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    /* sort up to position u + 1 where u = min(k, n_cols - 1). */
    size_t i;
    size_t sort_up_to = n_cols - 1;
    if (m_rows < n_cols) {
        sort_up_to = m_rows;
    }
    double sum_permanent = 0.0;
    double prod_permanent = 1.0;

    /* Allocate falling_fact, perm, and inv_perm arrays. */
    size_t falling_fact[128];
    size_t perm[128];
    size_t inv_perm[128];
    /* Set sentinal value. */
    falling_fact[0] = 0;

    /* Initialize the set to permute. */
    init_perm(n_cols, falling_fact, perm, inv_perm);
    // bool gen_next_perm();

    /* Handle first permutation. */
    for (i = 0; i < m_rows; ++i) {
            prod_permanent *= (ptr[i * n_cols + perm[i]]);
    }
    sum_permanent = prod_permanent;

    /* Iterate over second to last permutations. */
    while (gen_next_perm(falling_fact, perm, inv_perm, n_cols, sort_up_to)) {
        prod_permanent = 1.0;
        for (i = 0; i < m_rows; i++) {
            prod_permanent *= (ptr[i * n_cols + perm[i]]);
        }
        sum_permanent += prod_permanent;
    }
    return sum_permanent;
}


double glynn(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    size_t i, j;
    size_t pos = 0;
    size_t bound = n_cols - 1;
    size_t gray[64 + 1];

    int sign = 1;
    int delta[64];

    double sum, prod;
    double result = 1.0;
    double vec[64];

    /* Fill delta array (all +1 to start), and Gray code array with [0...n_cols]. */

    for (i = 0; i < n_cols; ++i) {
        gray[i] = i;
        delta[i] = 1;
    }

    /* Handle first Gray code. */

    for (j = 0; j < n_cols; ++j) {
        sum = 0.0;
        for (i = 0; i < n_cols; ++i) {
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

        /* Compute term. */

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


double glynn_rectangle(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    size_t i, j, k;
    size_t pos = 0;
    size_t bound = n_cols - 1;
    size_t gray[64 + 1];

    int sign = 1;
    int delta[64];

    double sum, prod;
    double result = 1.0;
    double vec[64];

    /* Fill delta array (all +1 to start), and Gray code array with [0...n_cols]. */

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

    /* Iterate from the second to the final Gray code. */

    while (pos != bound) {

        /* Update sign and delta. */

        sign *= -1;
        delta[bound - pos] *= -1;

        /* Compute term. */

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


double ryser(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    size_t i, j;
    unsigned long k;
    double rowsum, rowsumprod;
    double result = 0.0;

    /* Iterate over c = pow(2, n) submatrices (equal to (1 << n)) submatrices. */

    unsigned long c = 1UL << n_cols;

    /* Loop over columns of submatrix; compute product of row sums. */

    for (k = 0; k < c; ++k) {
        rowsumprod = 1.0;
        for (i = 0; i < n_cols; ++i) {

            /* Loop over rows of submatrix; compute row sum. */

            rowsum = 0.0;
            for (j = 0; j < n_cols; ++j) {

                /* Add element to row sum if the row index is in the characteristic *
                 * vector of the submatrix, which is the binary vector given by k.  */

                if (k & (1UL << j)) {
                    rowsum += ptr[n_cols * i + j];
                }
            }

            /* Update product of row sums. */

            rowsumprod *= rowsum;
        }

        /* Add term multiplied by the parity of the characteristic vector. */

        result += rowsumprod * (1 - ((__builtin_popcount(k) & 1) << 1));
    }

    /* Return answer with the correct sign (times -1 for odd n). */

    return result * ((n_cols % 2 == 1) ? -1 : 1);
}


double ryser_rectangle(const size_t m_rows, const size_t n_cols, double *const ptr)
{
    //return glynn_rectangle(m_rows, n_cols, ptr);
    /* Initialize all relevant variables. See combinatorial algorithm for more
     * details as it was already went over. */
        /* Initialize all relevant variables. See combinatorial algorithm for more
     * details as it was already went over. */
    size_t falling_fact[128];
    size_t perm[128];
    size_t inv_perm[128];
    falling_fact[0] = 0;

    size_t i, j, k, counter, sort_up_to;
    
    int value_sign = 1;
    
    double bin_c, sum_of_matrix_vals, prod_of_cols, result;
    double sum_over_k_vals = 0.0;
    double vec[64];
    
    /* Dealing with a rectangle. Can't use bit hacking trick here. */
    for (k = 0; k < m_rows; ++k) {
        /* Store the binomial coefficient for this k value bin_c. */
        bin_c = (double)(BINOMIAL((n_cols - m_rows + k), k));
        counter = 0; // Count how many permutations you have generated
        sum_of_matrix_vals = 0.0;
        prod_of_cols = 1.0;
        result = 0.0;
    
        /* (Re)initialize the set to permute for this k value. */
        init_perm(n_cols, falling_fact, perm, inv_perm);
    
        /* sort up to position u + 1 where u = min(m_rows - k, n_cols - 1). */
        sort_up_to = n_cols - 1;
    
        if ((m_rows - k) < sort_up_to) {
            sort_up_to = (m_rows - k);
        }
    
        for (i = 0; i < m_rows; ++i) {
            sum_of_matrix_vals = 0.0;
            for (j = 0; j < sort_up_to; ++j) {
                sum_of_matrix_vals += ptr[i * n_cols + perm[j]];
            }
            vec[i] = sum_of_matrix_vals;
            sum_of_matrix_vals = 0.0;
        }
        for (i = 0; i < m_rows; ++i) {
            prod_of_cols *= vec[i];
        }
    
        result += value_sign * bin_c * prod_of_cols;
        counter += 1;
    
        /* Iterate over second to last permutations of the set. */
        while (gen_next_perm(falling_fact, perm, inv_perm, n_cols, sort_up_to)) {
            sum_of_matrix_vals = 0.0;
            for (i = 0; i < m_rows; ++i) {
                sum_of_matrix_vals = 0.0;
                for (j = 0; j < m_rows - k; ++j) {
                    sum_of_matrix_vals += ptr[i * n_cols + perm[j]];
                }
                vec[i] = sum_of_matrix_vals;
                sum_of_matrix_vals = 0.0;
            }
            prod_of_cols = 1.0;
            for (i = 0; i < m_rows; ++i) {
                prod_of_cols *= vec[i];
            }
    
            result += value_sign * bin_c * prod_of_cols;
            counter += 1;
        }
        sum_over_k_vals += (result / (double)counter) * BINOMIAL(n_cols, (m_rows - k));
        value_sign *= -1;
    }
    return sum_over_k_vals;
}
