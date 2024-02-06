#include <stdio.h>
#include <stdlib.h>

#include "perm-mv0.h"
#include "kperm-gray.h"

#include "tables.h"

#include "permanent.h"


#ifdef WITH_TUNING_FILE

/* Include tuning file. */

#include "tuning.h"

#else

/* Set default tuning parameters. */

#define PARAM_1 1.0
#define PARAM_2 1.0
#define PARAM_3 1.0

#endif


double opt(const size_t m, const size_t n, double *const ptr)
{
    /* Use the fastest algorithm. */

    if (m < PARAM_1) {
        return (m == n) ? combinatoric(m, n, ptr) : combinatoric_rectangular(m, n, ptr);
    } else if (m * n < PARAM_2) {
        return (m == n) ? glynn(m, n, ptr) : glynn_rectangular(m, n, ptr);
    }
    else {
        return (m == n) ? ryser(m, n, ptr) : ryser_rectangular(m, n, ptr);
    }
}


double combinatoric(const size_t m, const size_t n, double *const ptr)
{
    (void)n;

    perm_mv0 permutations(m);

    const size_t *perm = permutations.data();

    size_t i;
    double prod;
    double out = 0.0;

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


double combinatoric_rectangular(const size_t m, const size_t n, double *const ptr)
{
    kperm_gray permutations(n);
    permutations.first(m);

    const size_t *perm = permutations.data();

    size_t i;
    double prod;
    double out = 0.0;

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


double glynn(const size_t m, const size_t n, double *const ptr)
{
    (void)n;

    size_t i, j;
    size_t pos = 0;
    size_t bound = m - 1;
    size_t perm[64 + 1];

    int sign = 1;
    int delta[64];

    double sum, prod;
    double out = 1.0;
    double vec[64];

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


double glynn_rectangular(const size_t m, const size_t n, double *const ptr)
{
    size_t i, j, k;
    size_t pos = 0;
    size_t bound = n - 1;
    size_t perm[64 + 1];

    int sign = 1;
    int delta[64];

    double sum, prod;
    double out = 1.0;
    double vec[64];

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



double ryser(const size_t m, const size_t n, double *const ptr)
{
    (void)n;

    size_t i, j;
    size_t k;
    double rowsum, rowsumprod;
    double out = 0;

    /* Iterate over c = pow(2, m) submatrices (equal to (1 << m)) submatrices. */

    size_t c = 1UL << m;

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


double ryser_rectangular(const size_t m, const size_t n, double *const ptr)
{
    kperm_gray permutations(n);

    const size_t *perm = permutations.data();

    size_t i, j, k, bin;
    int sign = 1;

    double colprod, matsum, permsum;
    double out = 0.0;
    double vec[64];

    /* Iterate over subsets from size 0 to size m */

    for (k = 0; k < m; ++k) {

        permutations.first(m - k);

        bin = BINOMIAL((n - m + k), k);
        permsum = 0.0;

        do {

            /* Compute permanents of each submatrix */

            for (i = 0; i < m; ++i) {
                matsum = 0.0;
                for (j = 0; j < m - k; ++j)
                    matsum += ptr[i * n + perm[j]];
                vec[i] = matsum;
            }

            colprod = 1.0;
            for (i = 0; i < m; ++i) {
                colprod *= vec[i];
            }

            permsum += colprod;
        }
        while (permutations.next());

        /* Add term to result */

        out += permsum * sign * bin;

        sign *= -1;
    }

    return out;
}
