/* The permanent commonly appears in problems related to quantum mechanics, and the most common
brute-force combinatorical method has time complexity O(N!N), thus it is useful to look
for more efficient algorithms. The two algorithms considered to be the fastest are one by
Ryser (based on the inclusion-exclusion principle), and one by Glynn (based on invariant theroy). 
All algorithms work for square NxN matrices, and are generalizable for MxN matrices.

The goal is to optimize the code and find the best algorithm for each value of M and N, and 
have a C++ function that will automatically find the best algorithm based on the size of the 
input matrix.

This code works with an input matrix in the form of Python NumPy array.

The function for generating k-permtations lexicographically was adapted from
"Matters Computational; Ideas, Algorithms, Source Code" by Jorg Arndt" Ch 12.1*/
/* ******** */


/* C headers. */
#include <stdbool.h>
#include <stdio.h>
#include <tgmath.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>


/* C functions. */
/* ********************************* */


/* A function to generate a table of binomial coefficients to be used in the Ryser
formula for rectangular matrices. Maximum values of N and K are set to 15 since anything larger is unnecessary
for our purposes. Recursion inspired by Pascal's triangle. We pass in the array as a constant
pointer and change the values inside the array. */

void bin_coeff(const uint64_t N, const uint64_t K, int64_t C[const N][K])
{
    for (int64_t k = 1; k <= 15; k++)
    {
        C[0][k] = 0;
    }
    for (int64_t n = 0; n <= 15; n++)
    {
        C[n][0] = 1;
    }

    for (int64_t n = 1; n <= 15; n++)
    {
        for (int64_t k = 1; k <= 15; k++)
        {
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
        }
    }
}


/* A function for swapping the values of two entries in an array while maintaining value of the pointer. */
void swap2(int64_t *perm, int64_t i, int64_t j)
{
    const int64_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}


/* A function that will initialize the set to permute. This first set is used as the first permutation. The
set is generated in ascedning order to ensure there are no smaller permutations possible. We keep track of the 
inverse permutation in order to simplify the swap update when generating the next permutation. The ffactorial
set is initialized to all zeroes and is used to know when we have generated all possible permutations.*/
void init_perm(const int64_t N, int64_t *const the_fact_set, int64_t *const the_perm_set, int64_t *const the_inv_set)
{
    for (int64_t i = 0; i < N; i++)
    {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}

/* A function that will use the current state of the permutation perm_ and update it to reflect that of the next permutation.
If the largest permutation has been generated, the function will return false, else it will determine the next permutation and 
update all three arrays (falling_fact, perm_, inv_perm). This was adapted from Jorgs "Matters Computational" C++ reference.*/
bool gen_next_perm(int64_t *const falling_fact, int64_t *const perm_, int64_t *const inv_perm, const int64_t rows, const int64_t cols, const int64_t u_)
{
    /* Use the current permutation to check for next permutation lexicographically, if it exists
    update the curr_array by swapping the leftmost changed element (at position i < k).
    Replace the elements up to position k by the smallest element lying to the right of i.
    Do not put remaining k, .., n - 1 positions in ascending order to improve efficiency
    has time complexity O(n) in the worst case. */
    int64_t i = rows - 1;
    int64_t m1 = cols - i - 1;
    /* begin update of falling_fact - recall falling_fact[0] = 0, so increment index accordingly.
    If i becomes negative during the check, you are pointing to the sentinal value, so you are
    done generating permutations. */
    while (falling_fact[i + 1] == m1)
    {
        falling_fact[i + 1] = 0;
        ++m1;
        --i;
    }
    if (i == -1)
    {
        return false;
    }
    ++falling_fact[i + 1];
    /* Find smallest element perm_[i] < perm_[j] that lies to the right of pos i,
    and then update the state of the permutation using its inverse to generate next. */
    int64_t z = (int64_t)perm_[i];
    do 
    {
        ++z;
    } while (inv_perm[z] <= i);
    const int64_t j = inv_perm[z];
    swap2(perm_, i, j);
    swap2(inv_perm, perm_[i], perm_[j]);
    ++i;
    int64_t y = 0;
    /* Find the smallest elements to the right of position i. */
    while (i < u_)
    {
        while (inv_perm[y] < i)
        {
            ++y;
        }
        const int64_t j = inv_perm[y];
        swap2(perm_, i, j);
        swap2(inv_perm, perm_[i], perm_[j]);
        ++i;
    }
    return true;
}


double combinatoric(double *const matrix, const int64_t m, const int64_t n)
{
    /* Store the pointer to the array in ptr. */
    double *ptr = matrix;

    /* Return the permanent of the matrix. */
    int64_t m_rows = m;
    int64_t n_cols = n;
    int64_t sort_up_to = n_cols - 1;
    /* sort up to position u + 1 where u = min(k, n_cols - 1). */
    if (m_rows < sort_up_to)
    {
        sort_up_to = m_rows;
    }
    double sum_permanent = 0.0;
    double prod_permanent = 1.0;
    /* Allocate falling_fact, perm_, and inv_perm arrays. */
    int64_t falling_fact[128];
    /* Set sentinal value. */
    falling_fact[0] = 0;
    int64_t perm_[128];
    int64_t inv_perm[128];
    init_perm(n_cols, falling_fact, perm_, inv_perm);
    bool gen_next_perm();
    /* Handle first permutation. */
    for (int64_t i = 0; i < m_rows; i++)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
    sum_permanent = prod_permanent;
    /* Iterate over second to last permutations. */
    while (gen_next_perm(falling_fact, perm_, inv_perm, m_rows, n_cols, sort_up_to))
    {
        prod_permanent = 1.0;
        for (int64_t i = 0; i < m_rows; i++)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
        sum_permanent += prod_permanent;
    }
    return sum_permanent;
}

/* This code for generating the permanent of a square NxN matrix was adapted + error corrected from Michael Richer's sketch code. */
double glynn(double *const matrix, const int64_t m, const int64_t n)
{
    /* Store the pointer to the array in ptr. */
    double *ptr = matrix;

    /* Return the permanent of the matrix. */
    int64_t m_rows = m;
    int64_t n_cols = n;

    /* Initialize gray code. */
    int64_t pos = 0;
    int64_t sign = 1;
    int64_t bound = n_cols - 1;
    int64_t delta[128];
    int64_t gray[128];

    /* Allocate and fill delta array (all +1 to start), and gray array from 0 to n_cols. */
    for (int64_t i = 0; i < n_cols; i++)
    {
        delta[i] = 1;
        gray[i] = i;
    }

    /* Allocate matrix for result of manual multiplies. */
    double vec[128];

    /* Handle first Gray code. */
    double result = 0.0, sum_ = 0.0, prod = 1.0;
    for (int64_t i = 0; i < n_cols; i++)
    {
        for (int64_t j = 0; j < m_rows; j++)
        {
            sum_ += (ptr[j * n_cols + i] * (double)delta[i]);
        }
        vec[i] = sum_;
    }

    for (int64_t i = 0; i < n_cols; i++)
    {
        prod *= vec[i];
    }
    for (int64_t i = 0; i < n_cols; i++)
    {
        result += (double)delta[i] * prod;
    }

    /* Update position and Gray code. */
    *gray = 0;
    *(gray + pos) = *(gray + pos + 1);
    ++pos;
    *(gray + pos) = pos;
    pos = gray[0];

    /* Iterate over second to last Gray codes. */
    while (pos != bound) 
    {
        /* Update sign and delta. */
        sign *= -1;
        *(delta + bound - pos) *= sign;
        /* Compute each Gray code term. */
        sum_ = 0.0;
        for (int64_t i = 0; i < n_cols; i++)
        {
            for (int64_t j = 0; j < m_rows; j++)
            {
                sum_ += (ptr[j * n_cols + i] * (double)delta[i]);
            }
            vec[i] = sum_;
        }
        prod = 1.0;
        for (int64_t i = 0; i < n_cols; i++)
        {
            prod *= vec[i];
        }
        for (int64_t i = 0; i < n_cols; i++)
        {
            result += (double)delta[i] * prod;
        }
        /* Go to next Gray code. */
        *gray = 0;
        *(gray + pos) = *(gray + pos + 1);
        ++pos;
        *(gray + pos) = pos;
        pos = gray[0];
    }

    /* Divide by external factor and return permanent. */
    return result / pow(2.0, (double)bound);
}


double ryser(double *const matrix, const int64_t m, const int64_t n)
{
    /* Store the pointer to the array in ptr. */
    double *ptr = matrix;

    /* Generate the binomial coefficient table. */
    int64_t C[20][20];
    void bin_coeff();
    bin_coeff(20, 20, C);

    /* Return the permanent of the matrix. */
    int64_t m_rows = m;
    int64_t n_cols = n;
    double sum_over_k_vals = 0.0;
    int64_t falling_fact[128];
    falling_fact[0] = 0;
    int64_t perm_[128];
    int64_t inv_perm[128];
    double vec[128];

    for (int64_t k = 0; k <= m_rows - 1; k++)
    {
        /* Store the binomial coefficient for this k value bin_c. */
        int64_t bin_c = C[n_cols - m_rows + k][k];

        double sum_of_matrix_vals = 0.0;
        double prod_of_cols = 1.0;
        double result = 0.0;
        double value_sign = pow(-1.0, (double)k);
        /* (Re)initialize the set to permute for this k value. */
        init_perm(n_cols, falling_fact, perm_, inv_perm);
        bool gen_next_perm();
        /* sort up to position u + 1 where u = min(k, n_cols - 1). */
        int64_t sort_up_to = n_cols - 1;
        if (m_rows - k < sort_up_to)
        {
            sort_up_to = m_rows - k;
        }

        for (int64_t i = 0; i < m_rows; i++)
        {
            for (int64_t j = 0; j < m_rows - k; j++)
            {
                sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
            }
            vec[i] = sum_of_matrix_vals;
        }
        for (int64_t i = 0; i < m_rows; i++)
        {
            prod_of_cols *= vec[i];
        }
        result += value_sign * (double)bin_c * prod_of_cols;
        
        /* Iterate over second to last permutations of the set. */
        while (gen_next_perm(falling_fact, perm_, inv_perm, m_rows - k, n_cols, sort_up_to))
        {
            sum_of_matrix_vals = 0.0;
            for (int64_t i = 0; i < m_rows; i++)
            {
                for (int64_t j = 0; j < m_rows - k; j++)
                {
                    sum_of_matrix_vals += (ptr[i * n_cols + perm_[j]]);
                }
                vec[i] = sum_of_matrix_vals;
            }
            prod_of_cols = 1.0;
            for (int64_t i = 0; i < m_rows; i++)
            {
                prod_of_cols *= vec[i];
            }
            result += value_sign * (double)bin_c * prod_of_cols;
        }
        sum_over_k_vals += result;
    }
    return sum_over_k_vals;
}

enum winning_algorithm
{
    COMBINATORIC,
    RYSER,
    GLYNN
};

int main(void)
{
    /* Initialize maximum rows and columns of interest. */
    int64_t M = 3, N = 3;

    if (M > 15 || N > 15)
    {
        printf("Matrix dimensions must be less than 15! Larger values are too computationally intensive!! \n");
        printf("This program computes the permanent of a (not necessarily sqaure) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 15, since anything larger is too intensive. \nint64_t m = M, n = N ------ where M is the largest number of rows of interest, and N is the largest number of columns of interest. \nExample call: \nint64_t m = 2, n = 2; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \nand print out the fastest algorithm for that matrix size.\n");
        return -1;
    }
    else if (M > N)
    {
        printf("Program is designed to take matrix dimensions M <= N where M = # rows and N = # cols. Per(A) = Per(A_transpose), please transpose your matrix! \n");
        printf("This program computes the permanent of a (not necessarily sqaure) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 15, since anything larger is too intensive. \nint64_t m = M, n = N ------ where M is the largest number of rows of interest, and N is the largest number of columns of interest. \nExample call: \nint64_t m = 2, n = 2; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \nand print out the fastest algorithm for that matrix size.\n");
        return -1;
    }
    else if (!M || !N)
    {
        printf("You have entered an empty matrix! The permanent has no meaning. Come back when you need a permanent :) \n");
        printf("This program computes the permanent of a (not necessarily sqaure) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 15, since anything larger is too intensive. \nint64_t m = M, n = N ------ where M is the largest number of rows of interest, and N is the largest number of columns of interest. \nExample call: \nint64_t m = 2, n = 2; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \nand print out the fastest algorithm for that matrix size.\n");
        return -1;
    }
    else if (M < 0 || N < 0)
    {
        printf("You have given a negative value for a matrix dimension! Dimensions must be positive values. \n");
        printf("This program computes the permanent of a (not necessarily sqaure) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 15, since anything larger is too intensive. \nint64_t m = M, n = N ------ where M is the largest number of rows of interest, and N is the largest number of columns of interest. \nExample call: \nint64_t m = 2, n = 2; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \nand print out the fastest algorithm for that matrix size.\n");
        return -1;
    }
    else if (M == 1 && N == 1)
    {
        printf("The permanent of a 1x1 matrix/single value is the value itself! Enter something interesting to compute. \n");
        printf("This program computes the permanent of a (not necessarily sqaure) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 15, since anything larger is too intensive. \nint64_t m = M, n = N ------ where M is the largest number of rows of interest, and N is the largest number of columns of interest. \nExample call: \nint64_t m = 2, n = 2; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \nand print out the fastest algorithm for that matrix size.\n");
        return -1;
    }

    for (int64_t m = 2; m <= M; m++)
    {
        for (int64_t n = m; n <= N; n++)
        {
            /* Find the winner 10 times to account for background programs slowing down computation. */
            int64_t combinatoric_counter = 0;
            int64_t glynn_counter = 0;
            int64_t ryser_counter = 0;

            for (int64_t i = 0; i < 10; i++)
            {

                double matrix[256];
                double counter = 1.0;
                for (int64_t i = 0; i < m; i++)
                {
                    for (int64_t j = 0; j < n; j++)
                    {
                        matrix[i * n + j] = counter;
                        counter += 1.0;
                    }
                }
                clock_t begin_1 = clock();
                double solve_combinatoric = combinatoric((double *)matrix, m, n);
                clock_t end_1 = clock();
                double time_spent_on_comb = (double)(end_1 - begin_1) / CLOCKS_PER_SEC;
                printf("The permanent is: %f", solve_combinatoric);
                clock_t begin_2 = clock();
                double solve_glynn = glynn((double *)matrix, m, n);
                clock_t end_2 = clock();
                printf("The permanent is: %f", solve_glynn);
                double time_spent_on_glynn = (double)(end_2 - begin_2) / CLOCKS_PER_SEC;
                clock_t begin_3 = clock();
                double solve_ryser = ryser((double *)matrix, m, n);
                clock_t end_3 = clock();
                printf("The permanent is: %f", solve_ryser);
                double time_spent_on_ryser = (double)(end_3 - begin_3) / CLOCKS_PER_SEC;
                printf("\n Time spent on combinatoric algorithm: %f. Time spent on Glynn algorithm: %f. \
                Time spent on Ryser algorithm: %f.", time_spent_on_comb, time_spent_on_glynn, time_spent_on_ryser);

                if (solve_combinatoric <= solve_ryser && solve_combinatoric <= solve_glynn)
                {
                    combinatoric_counter += 1;
                }
                else if (solve_glynn <= solve_ryser)
                {
                    glynn_counter += 1;
                }
                else 
                {
                    ryser_counter += 1;
                }
            }
            enum winning_algorithm alg;
            if (combinatoric_counter > glynn_counter && combinatoric_counter > ryser_counter)
            {
                alg = COMBINATORIC; 
            }
            else if (glynn_counter > ryser_counter)
            {
                alg = GLYNN;
            }
            else
            {
                alg = RYSER;
            }
            switch (alg)
            {
                case COMBINATORIC:
                    printf("For m/n=%1.4f, m+n=%lld, the fastest algorithm is the combinatoric algorithm. \n", (double)m/(double)n, m+n);
                    break;
                case GLYNN:
                    printf("For m/n=%1.4f, m+n=%lld, the fastest algorithm is the glynn algorithm. \n", (double)m/(double)n, m+n);
                    break;
                case RYSER:
                    printf("For m/n=%1.4f, m+n=%lld, the fastest algorithm is the ryser algorithm. \n", (double)m/(double)n, m+n);
                    break;
            }   
        }
    }
}
