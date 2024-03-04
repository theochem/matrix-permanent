#ifndef PERMANENT_PERMANENT_H
#define PERMANENT_PERMANENT_H


#include <cstdlib>

#include <complex>

#include "perm-mv0.h"
#include "kperm-gray.h"

#include "tables.h"


#ifdef WITH_TUNING_FILE

/* Include tuning file. */

#include "tuning.h"

#else

/* Set default tuning parameters. */

#define PARAM_1 -0.572098
#define PARAM_2 -22.014212
#define PARAM_3 15.297940
#define PARAM_4 3.0

#endif /* WITH_TUNING_FILE */


/* Allow type promotion to complex. */

template <typename T>
struct identity_t { typedef T type; };

#define COMPLEX_OPS(OP) \
    \
    template <typename _Tp> \
    std::complex<_Tp> \
    operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
    { return lhs OP rhs; } \
    \
    template <typename _Tp> \
    std::complex<_Tp> \
    operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
    { return lhs OP rhs; } \

COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)

#undef COMPLEX_OPS


template<typename T>
T combinatoric(const std::size_t m, const std::size_t n, T *const ptr)
{
    (void)n;

    perm_mv0 permutations(m);

    const std::size_t *perm = permutations.data();

    std::size_t i;
    T prod;
    T out = 0.0;

    do {

        prod = 1.0;

        for (i = 0; i < m; ++i) {
            prod *= ptr[i * m + perm[i]];
        }

        out += prod;

    }
    while (permutations.next());

    return out;
}


template<typename T>
T combinatoric_rectangular(const std::size_t m, const std::size_t n, T *const ptr)
{
    kperm_gray permutations(n);
    permutations.first(m);

    const std::size_t *perm = permutations.data();

    std::size_t i;
    T prod;
    T out = 0.0;

    do {

        prod = 1.0;

        for (i = 0; i < m; ++i) {
            prod *= ptr[i * n + perm[i]];
        }

        out += prod;

    }
    while (permutations.next());

    return out;
}


template<typename T>
T glynn(const std::size_t m, const std::size_t n, T *const ptr)
{
    (void)n;

    std::size_t i, j;
    std::size_t pos = 0;
    std::size_t bound = m - 1;
    std::size_t perm[64 + 1];

    int sign = 1;
    int delta[64];

    T sum, prod;
    T out = 1.0;
    T vec[64];

    /* Fill delta array (all +1 to start), and permutation array with [0...m]. */

    for (i = 0; i < m; ++i) {
        perm[i] = i;
        delta[i] = 1;
    }

    /* Handle first permutation. */

    for (j = 0; j < m; ++j) {
        sum = 0.0;
        for (i = 0; i < m; ++i) {
            sum += ptr[i * m + j] * delta[i];
        }
        vec[j] = sum;
    }
    for (j = 0; j < m; ++j) {
        out *= vec[j];
    }

    /* Iterate from the second to the final permutation. */

    while (pos != bound) {

        /* Update sign and delta. */

        sign *= -1;
        delta[bound - pos] *= -1;

        /* Compute term. */

        for (j = 0; j < m; ++j) {
            sum = 0.0;
            for (i = 0; i < m; ++i) {
                sum += ptr[i * m + j] * delta[i];
            }
            vec[j] = sum;
        }

        /* Multiply by the product of the vectors in delta. */

        prod = 1.0;
        for (i = 0; i < m; ++i) {
            prod *= vec[i];
        }
        out += sign * prod;

        /* Go to next permutation. */

        perm[0] = 0;
        perm[pos] = perm[pos + 1];
        ++pos;
        perm[pos] = pos;
        pos = perm[0];
    }

    /* Divide by external factor and return permanent. */

    return out / (1UL << bound);
}


template<typename T>
T glynn_rectangular(const std::size_t m, const std::size_t n, T *const ptr)
{
    std::size_t i, j, k;
    std::size_t pos = 0;
    std::size_t bound = n - 1;
    std::size_t perm[64 + 1];

    int sign = 1;
    int delta[64];

    T sum, prod;
    T out = 1.0;
    T vec[64];

    /* Fill delta array (all +1 to start), and permutation array with [0...n]. */

    for (i = 0; i < n; ++i) {
        delta[i] = 1;
        perm[i] = i;
    }

    /* Handle first permutation. */

    for (j = 0; j < n; ++j) {
        sum = 0.0;
        for (i = 0; i < m; ++i) {
            sum += ptr[i * n + j] * delta[i];
        }
        for (k = m; k < n; ++k) {
            sum += delta[k];
        }
        vec[j] = sum;
    }
    for (i = 0; i < n; ++i) {
        out *= vec[i];
    }

    /* Iterate from the second to the final permutation. */

    while (pos != bound) {

        /* Update sign and delta. */

        sign *= -1;
        delta[bound - pos] *= -1;

        /* Compute term. */

        for (j = 0; j < n; ++j) {
            sum = 0.0;
            for (i = 0; i < m; ++i) {
                sum += ptr[i * n + j] * delta[i];
            }

            for (k = m; k < n; ++k) {
                sum += delta[k];
            }
            vec[j] = sum;
        }

        /* Multiply by the product of the vectors in delta. */

        prod = 1.0;
        for (i = 0; i < n; ++i) {
            prod *= vec[i];
        }
        out += sign * prod;

        /* Go to next permutation. */

        perm[0] = 0;
        perm[pos] = perm[pos + 1];
        ++pos;
        perm[pos] = pos;
        pos = perm[0];
    }

    /* Divide by external factor and return permanent. */

    return out / ((1UL << bound) * FACTORIAL(n - m));
}


template<typename T>
T ryser(const std::size_t m, const std::size_t n, T *const ptr)
{
    (void)n;

    std::size_t i, j;
    std::size_t k;
    T rowsum, rowsumprod;
    T out = 0;

    /* Iterate over c = pow(2, m) submatrices (equal to (1 << m)) submatrices. */

    std::size_t c = 1UL << m;

    /* Loop over columns of submatrix; compute product of row sums. */

    for (k = 0; k < c; ++k) {
        rowsumprod = 1.0;
        for (i = 0; i < m; ++i) {

            /* Loop over rows of submatrix; compute row sum. */

            rowsum = 0.0;
            for (j = 0; j < m; ++j) {

                /* Add element to row sum if the row index is in the characteristic *
                 * vector of the submatrix, which is the binary vector given by k.  */

                if (k & (1UL << j)) {
                    rowsum += ptr[m * i + j];
                }
            }

            /* Update product of row sums. */

            rowsumprod *= rowsum;
        }

        /* Add term multiplied by the parity of the characteristic vector. */

        out += rowsumprod * (1 - ((__builtin_popcountll(k) & 1) << 1));
    }

    /* Return answer with the correct sign (times -1 for odd m). */

    return out * ((m % 2 == 1) ? -1 : 1);
}


// template<typename T>
// T ryser_rectangular(const std::size_t m, const std::size_t n, T *const ptr)
// {
//     kperm_gray permutations(n);
//
//     const std::size_t *perm = permutations.data();
//
//     std::size_t i, j, k, bin;
//     int sign = 1;
//
//     T colprod, matsum, permsum;
//     T out = 0.0;
//     T vec[64];
//
//     /* Iterate over subsets from size 0 to size m */
//
//     for (k = 0; k < m; ++k) {
//
//         permutations.first(m - k);
//
//         bin = BINOMIAL((n - m + k), k);
//         permsum = 0.0;
//
//         do {
//
//             /* Compute permanents of each submatrix */
//
//             for (i = 0; i < m; ++i) {
//                 matsum = 0.0;
//                 for (j = 0; j < m - k; ++j)
//                     matsum += ptr[i * n + perm[j]];
//                 vec[i] = matsum;
//             }
//
//             colprod = 1.0;
//             for (i = 0; i < m; ++i) {
//                 colprod *= vec[i];
//             }
//
//             permsum += colprod;
//         }
//         while (permutations.next());
//
//         /* Add term to result */
//
//         out += permsum * sign * bin;
//
//         sign *= -1;
//     }
//
//     return out;
// }


namespace {

inline void swap2(std::size_t *perm, std::size_t i, std::size_t j)
{
    const std::size_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}


void init_perm(const std::size_t N, std::size_t *const the_fact_set, std::size_t *const the_perm_set, std::size_t *const the_inv_set)
{
    for (std::size_t i = 0; i < N; i++) {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}


bool gen_next_perm(std::size_t *const falling_fact, std::size_t *const perm, std::size_t *const inv_perm, const std::size_t cols, const std::size_t u_)
{
    /* Use the current permutation to check for next permutation
     * lexicographically, if it exists update the curr_array by swapping the
     * leftmost changed element (at position i < k). Replace the elements up to
     * position k by the smallest element lying to the right of i. Do not put
     * remaining k, .., n - 1 positions in ascending order to improve efficiency
     * has time complexity O(n) in the worst case. */
    std::size_t i = u_ - 1;
    std::size_t m1 = cols - i - 1;

    /* Begin update of falling_fact - recall falling_fact[0] = 0, so increment
     * index accordingly. If i becomes negative during the check, you are
     * pointing to the sentinel value, so you are done generating
     * permutations. */
    while (falling_fact[i + 1] == m1) {
        falling_fact[i + 1] = 0;
        ++m1;
        --i;
    }
    if (i == 0UL - 1) {
        return false;
    }
    ++falling_fact[i + 1];

    /* Find smallest element perm[i] < perm[j] that lies to the right of
     * pos i, and then update the state of the permutation using its inverse
     * to generate next. */
    std::size_t z = (std::size_t)perm[i];
    do {
        ++z;
    } while (inv_perm[z] <= i);
    const std::size_t j = inv_perm[z];
    swap2(perm, i, j);
    swap2(inv_perm, perm[i], perm[j]);
    ++i;
    std::size_t y = 0;

    /* Find the smallest elements to the right of position i. */
    while (i < u_) {
        while (inv_perm[y] < i) {
            ++y;
        }
        const std::size_t j = inv_perm[y];
        swap2(perm, i, j);
        swap2(inv_perm, perm[i], perm[j]);
        ++i;
    }
    return true;
}

} /* namespace */


template<typename T>
T ryser_rectangular(const std::size_t m, const std::size_t n, T *const ptr)
{
    std::size_t falling_fact[128];
    std::size_t perm[128];
    std::size_t inv_perm[128];
    falling_fact[0] = 0;

    std::size_t i, j, k, counter, sort_up_to;

    int value_sign = 1;

    T bin_c, sum_of_matrix_vals, prod_of_cols, result;
    T sum_over_k_vals = 0.0;
    T vec[64];

    /* Dealing with a rectangle. Can't use bit hacking trick here. */
    for (k = 0; k < m; ++k) {
        /* Store the binomial coefficient for this k value bin_c. */
        bin_c = (double)(BINOMIAL((n - m + k), k));
        counter = 0; // Count how many permutations you have generated
        sum_of_matrix_vals = 0.0;
        prod_of_cols = 1.0;
        result = 0.0;

        /* (Re)initialize the set to permute for this k value. */
        init_perm(n, falling_fact, perm, inv_perm);

        /* sort up to position u + 1 where u = min(m - k, n - 1). */
        sort_up_to = n - 1;

        if ((m - k) < sort_up_to) {
            sort_up_to = (m - k);
        }

        for (i = 0; i < m; ++i) {
            sum_of_matrix_vals = 0.0;
            for (j = 0; j < sort_up_to; ++j) {
                sum_of_matrix_vals += ptr[i * n + perm[j]];
            }
            vec[i] = sum_of_matrix_vals;
            sum_of_matrix_vals = 0.0;
        }
        for (i = 0; i < m; ++i) {
            prod_of_cols *= vec[i];
        }

        result += value_sign * bin_c * prod_of_cols;
        counter += 1;

        /* Iterate over second to last permutations of the set. */
        while (gen_next_perm(falling_fact, perm, inv_perm, n, sort_up_to)) {
            sum_of_matrix_vals = 0.0;
            for (i = 0; i < m; ++i) {
                sum_of_matrix_vals = 0.0;
                for (j = 0; j < m - k; ++j) {
                    sum_of_matrix_vals += ptr[i * n + perm[j]];
                }
                vec[i] = sum_of_matrix_vals;
                sum_of_matrix_vals = 0.0;
            }
            prod_of_cols = 1.0;
            for (i = 0; i < m; ++i) {
                prod_of_cols *= vec[i];
            }

            result += value_sign * bin_c * prod_of_cols;
            counter += 1;
        }
        sum_over_k_vals += (result / (double)counter) * BINOMIAL(n, (m - k));
        value_sign *= -1;
    }
    return sum_over_k_vals;
}


template<typename T>
T opt(const std::size_t m, const std::size_t n, T *const ptr)
{
    /* Use the fastest algorithm. */

    if (m == n && n <= PARAM_4) {
        return ryser(m, n, ptr);
    } else if (n * PARAM_1 + ((double)m / n) * PARAM_2 + PARAM_3 > 0) {
        return (m == n) ? combinatoric(m, n, ptr) : combinatoric_rectangular(m, n, ptr);
    } else {
        return (m == n) ? glynn(m, n, ptr) : glynn_rectangular(m, n, ptr);
    }
}


#endif /* PERMANENT_PERMANENT_H */
