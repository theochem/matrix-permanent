/* ******************************************** */
/* The permanent commonly appears in problems related to quantum mechanics, and the most common brute-force combinatorial method has time complexity O(N!N), thus it is useful to look for more efficient algorithms. The two algorithms considered to be the fastest are one by Ryser (based on the inclusion-exclusion principle), and one by Glynn (based on invariant theory). All algorithms work for square NxN matrices, and are generalizable for MxN matrices.

The goal is to optimize the code and find the best algorithm for each value of M and N, and have a C++ function that will automatically find the best algorithm based on the size of the input matrix.

This code works by reading a user defined input file matrix_dimensions.csv of format M, N, r for reading specifying the matrix size with M rows and N columns respresenting the largest desired matrix dimensions for solving the permanent. The program will solve the permanent of all matrices of sizes mxn the specified number of times by the user where: m = 2; m <= M; m++; n = m; n <= N; n++; and write out to the file fast_permanent.csv data regarding fastest algorithm for that matrix size. Data corresponding to M/N, Size, Fastest Algorithm, M, N, Mean Time to Solve, Standard Deviation, xFaster Combinatorial, xFaster Glynn, xFaster Ryser, Speed Combinatorial, Speed Glynn, Speed Ryser will be written to the output file.

There is a provided Python code fastest_permanent_plotter.py that will visualize the data for your convenience. It will make it easier to decipher boundaries for which squareness and size combination each algorithm performs best.*/

/* ******************************************** */

/* C headers. */
/* ******************************************** */

#include <stdbool.h>
#include <stdio.h>
#include <tgmath.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>

/* C functions. */
/* ******************************************** */

/* A function to generate a table of binomial coefficients to be used in the Ryser formula for rectangular matrices. Maximum values of N and K are set to 10 since anything larger is unnecessary for our purposes. Recursion is inspired by Pascal's triangle. We pass in the array as a constant pointer and change the values inside the array. Many people have written this code. Mine turns out to be very similar to one found in a stackoverflow discussion: https://stackoverflow.com/questions/11032781/fastest-way-to-generate-binomial-coefficients */

void bin_coeff(const uint64_t N, const uint64_t K, int64_t C[const N][K])
{
    for (int64_t k = 1; k <= 10; k++)
    {
        /* Set the value to 0 when we are choosing from an empty set. */
        C[0][k] = 0;
    }
    for (int64_t n = 0; n <= 10; n++)
    {
        /* Set the value to 1 when we are choosing 0 items from a set. */
        C[n][0] = 1;
    }

    for (int64_t n = 1; n <= 10; n++)
    {
        for (int64_t k = 1; k <= 10; k++)
        {
            /* Use recursion and the pattern defined in Pascal's triangle to populate the table of binomial coefficients. */
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k];
        }
    }
}

/* A function for swapping the values of two entries in an array while maintaining value of the pointer. In class we went over a similar swapping function where the values stored in variables were swapped. Here matrix entries are swapped. */

void swap2(int64_t *perm, int64_t i, int64_t j)
{
    const int64_t temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
}

/* A function that will initialize the set to permute. This first set is used as the first permutation. The set is generated in ascending order to ensure there are no smaller permutations possible. We keep track of the inverse permutation in order to simplify the swap update when generating the next permutation. The ffactorial set is initialized to all zeroes and is used to know when we have generated all possible permutations. The idea for this function comes from the use of constructors and destructors in the class kperm_lex from "Matters Computational; Ideas, Algorithms, Source Code" by Jorg Arndt" Ch 12.1 https://www.jjj.de/fxt/fxtbook.pdf */

void init_perm(const int64_t N, int64_t *const the_fact_set, int64_t *const the_perm_set, int64_t *const the_inv_set)
{
    for (int64_t i = 0; i < N; i++)
    {
        the_fact_set[i + 1] = 0;
        the_perm_set[i] = i;
        the_inv_set[i] = i;
    }
}

/* A function that will use the current state of the permutation perm_ and update it to reflect that of the next permutation. If the largest permutation has been generated the function will return false, else it will determine the next permutation and update all three arrays (falling_fact, perm_, inv_perm). This function was adapted from the function next() in "Matters Computational; Ideas, Algorithms, Source Code" by Jorg Arndt" Ch 12.1 https://www.jjj.de/fxt/fxtbook.pdf */

bool gen_next_perm(int64_t *const falling_fact, int64_t *const perm_, int64_t *const inv_perm, const int64_t cols, const int64_t u_)
{
    /* Use the current permutation to check for next permutation lexicographically, if it exists update the curr_array by swapping the leftmost changed element (at position i < k). Replace the elements up to position k by the smallest element lying to the right of i. Do not put remaining k, .., n - 1 positions in ascending order to improve efficiency has time complexity O(n) in the worst case. */

    int64_t i = u_ - 1;
    int64_t m1 = cols - i - 1;

    /* Begin update of falling_fact - recall falling_fact[0] = 0, so increment index accordingly. If i becomes negative during the check, you are pointing to the sentinel value, so you are done generating permutations. */

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

    /* Find smallest element perm_[i] < perm_[j] that lies to the right of pos i, and then update the state of the permutation using its inverse to generate next. */

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

/* Solve the permanent using the combinatorial algorithm. */

double combinatorial(double *const matrix, const int64_t m, const int64_t n)
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
    falling_fact[0] = 0; // Set sentinel value.
    int64_t perm_[128];
    int64_t inv_perm[128];

    init_perm(n_cols, falling_fact, perm_, inv_perm); // Initialize the set to permute
    bool gen_next_perm();

    /* Handle first permutation. */

    for (int64_t i = 0; i < m_rows; i++)
        {
            prod_permanent *= (ptr[i * n_cols + perm_[i]]);
        }
    sum_permanent = prod_permanent;

    /* Iterate over second to last permutations. */

    while (gen_next_perm(falling_fact, perm_, inv_perm, n_cols, sort_up_to))
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

/* Solve the permanent using the Glynn algorithm. This trick used for updating the position and gray code was adapted + corrected from Michael Richer's sketch code (lines 262-293). He just sent it to me over Microsoft teams. */

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

    if (m_rows == n_cols) // Dealing with a square matrix
    {
        /* Handle first Gray code. */

        double result = 1.0;
        for (int64_t j = 0; j < n_cols; j++)
        {
            double sum_ = 0.0;
            for (int64_t i = 0; i < n_cols; i++)
            {
                /* Sum over all the values in each column. */
                sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
            }
            vec[j] = sum_;
        }

        for (int64_t j = 0; j < n_cols; j++)
        {
            result *= vec[j];
        }
        
        /* Iterate over second to last Gray codes. */

        while (pos != bound) 
        {
            /* Update sign and delta. */
            sign *= -1;
            *(delta + bound - pos) *= -1;

            /* Compute each Gray code term. */

            for (int64_t j = 0; j < n_cols; j++)
            {
                double sum_ = 0.0;
                for (int64_t i = 0; i < n_cols; i++)
                {
                    sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
                }
                vec[j] = sum_;
            }
            double prod = 1.0;
            for (int64_t i = 0; i < n_cols; i++)
            {
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

    else // Dealing with a rectangle
    {
        /* Handle first Gray code. */

        double result = 1.0;
        for (int64_t j = 0; j < n_cols; j++)
        {
            double sum_ = 0.0;
            for (int64_t i = 0; i < m_rows; i++)
            {
                sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
            }
            for (int64_t k = m_rows; k < n_cols; k++)
            {
                sum_ += (double)delta[k];
            }
            vec[j] = sum_;
        }

        for (int64_t i = 0; i < n_cols; i++)
        {
            result *= vec[i];
        }        

        /* Iterate over second to last Gray codes. */

        while (pos != bound) 
        {
            /* Update sign and delta. */
            sign *= -1;
            *(delta + bound - pos) *= -1;

            /* Compute each Gray code term. */

            for (int64_t j = 0; j < n_cols; j++)
            {
                double sum_ = 0.0;
                for (int64_t i = 0; i < m_rows; i++)
                {
                    sum_ += (ptr[i * n_cols + j] * (double)delta[i]);
                }

                for (int64_t k = m_rows; k < n_cols; k++)
                {
                    sum_ += (double)delta[k];
                }
                vec[j] = sum_;
            }
            double prod = 1.0;
            for (int64_t i = 0; i < n_cols; i++)
            {
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

/* Solve the permanent using the Ryser algorithm. */

double ryser(double *const matrix, const int64_t m, const int64_t n, int64_t *const binc)
{
    /* Store the pointer to the array in ptr. */

    double *ptr = matrix;

    /* Store the pointer to the binomial coefficient table in C. */

    int64_t *C = binc;

    /* Return the permanent of the matrix. */

    /* Initialize all relevant variables. See combinatorial algorithm for more details as it was already went over. */

    int64_t m_rows = m;
    int64_t n_cols = n;
    double sum_over_k_vals = 0.0;
    int64_t falling_fact[128];
    falling_fact[0] = 0;
    int64_t perm_[128];
    int64_t inv_perm[128];
    double vec[128];


    if (m_rows == n_cols) // Dealing with a square matrix. This bit-hacking trick was modified from C++ code from Micahel Richer (lines 393-428)
    {
        int32_t i, j, k;
        int64_t sum = 0, rowsum, rowsumprod;

        /* Iterate over c = pow(2, n) submatrices (equal to (1 << n)) submatrices. */ 
        int32_t c = 1UL << n_cols;

        /* Loop over columns of submatrix; compute product of row sums. */
        for (k = 0; k < c; k++)
        {
            rowsumprod = 1;
            for (i = 0; i < n_cols; i++)
            {
                /* Loop over rows of submatrix; compute row sum. */
                rowsum = 0;
                for (j = 0; j < n_cols; j++)
                {
                    /* Add element to row sum if the row index is in the characteristic vector of the submatrix, which is the binary vector given by k. */
                    if (k & (1UL << j))
                    {
                        rowsum += matrix[n_cols * i + j];
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
        if (n_cols % 2 == 1)
        {
            /* n is odd. */
            sign *= -1;
        }
        return sum * sign;
    }
    
    else // Dealing with a rectangle. Can't use bit hacking trick here.
    {

        int32_t value_sign = 1;
        for (int64_t k = 0; k < m_rows; k++)
        {
            /* Store the binomial coefficient for this k value bin_c. */

            int64_t bin_c = C[20 * (n_cols - m_rows + k) + k];

            double sum_of_matrix_vals = 0.0;
            double prod_of_cols = 1.0;
            double result = 0.0;

            /* (Re)initialize the set to permute for this k value. */

            init_perm(n_cols, falling_fact, perm_, inv_perm);
            bool gen_next_perm();

            /* sort up to position u + 1 where u = min(m_rows - k, n_cols - 1). */

            int64_t sort_up_to = n_cols - 1;

            if ((m_rows - k) < sort_up_to)
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

            while (gen_next_perm(falling_fact, perm_, inv_perm, n_cols, sort_up_to))
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
            value_sign *= -1;
        }
        
        return sum_over_k_vals;   
    }
}

enum winning_algorithm
{
    COMBINATORIAL,
    RYSER,
    GLYNN
};

int main(void)
{
    /* Use fscanf to take a file pointer to the file containing required information about maximum matrix dimensions and number of calculations to run, and read the data from the files in the required particular format. A user will have a csv file called "matrix_dimensions.csv" storing integers in the same location as this c file. This structure for opening and reading a file is adapted from lecture. */

    const char filename[] = "matrix_dimensions.csv";
    int num_read = 0;
    FILE *file_ptr = NULL;
    int64_t m = 0, n = 0, r = 0;
    uint64_t line = 1;

    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
        return -1;
    }

    char correct_input[] = "This program computes the permanent of a (not necessarily square) MxN matrix, where M <= N (since per(A) = per(A_transpose)). \nThe maximum values of M, N = 10, since anything larger is too computationally intensive. \nInput data types are as follows: \nint64_t M, N, r ------ where M is the largest number of rows of interest, N is the largest number of columns of interest, and r is the number of times to run each calculation. You must give a value of r between 10 and 100. \nExample csv: \n4, 6, 10; \nThe program will solve the permanent of all matrices of sizes mxn where \nm = 2; m <= M; m++; \n    n = m; n <= N; n++; \n r times for each algorithm and write information regarding computation time to a fast_permanent.csv file.";

    do
    {
        num_read = fscanf(file_ptr, "%" PRId64 ",%" PRId64 ",%" PRId64, &m, &n, &r);
        if (num_read == 3 && m == 1 && n == 1)
        {
            printf("The permanent of a 1x1 matrix/single value is the value itself! Enter something interesting to compute. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && m < 0 || n < 0)
        {
            printf("You have given a negative value for a matrix dimension! Dimensions must be positive values. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && !n || !m)
        {
            printf("You have entered an empty matrix! The permanent has no meaning. Come back when you need a permanent :) \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && m > n)
        {
            printf("Program is designed to take matrix dimensions M <= N where M = # rows and N = # cols. Per(A) = Per(A_transpose), please transpose your matrix! \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && m > 10 || n > 10)
        {
            printf("Matrix dimensions must be less than 10! Larger values are too computationally intensive!! \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && !r)
        {
            printf("Please specify how many times you would like to run each calculation. It is suggested that this number is at least 10 for statistical significance. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && r < 0)
        {
            printf("Please specify how many times you would like to run each calculation. This number must be at least 10 for statistical significance. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && r < 10)
        {
            printf("Please specify how many times you would like to run each calculation. This number must be at least 10 for statistical significance. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3 && r > 100)
        {
            printf("WARNING! This program does not allow you to run each calculation more than 100 times as it is meant to be a quick test to see which algorithms perform best for different values of M and N. Please give a value between 10 and 100. \n ------------------------------ \n");
            printf("%s\n", correct_input);
            return -1;
        }
        else if (num_read == 3)
        {
            printf("Read values: %" PRId64 ", %" PRId64 ", %" PRId64 " on line %" PRId64 ". Maximum matrix size will be %" PRId64 "-by-%" PRId64 ". \nThe permutations will be solved %" PRId64 " times for each algorithm.\n", m, n, r, line, m, n, r);
        }
        else
        {
            printf("Error: Encountered invalid data on line %" PRId64 "!\n", line);
            printf("%s\n", correct_input);
            return -1;
        }
        line++;
    } while (!feof(file_ptr));

    if (feof(file_ptr))
        printf("Successfully read the file %s (%" PRId64 " line total).\n", filename, line - 1);

    fclose(file_ptr);

    /* Open a file for writing to. This code was adapted from the structure we went over in lecture. */

    const char filename2[] = "fast_permanent.csv";
    FILE *file_ptr2 = NULL;

    file_ptr2 = fopen(filename2, "w");
    if (file_ptr2 == NULL)
    {
        perror("Oh no! Cannot open file!");
        return -1;
    }
    
    if (fprintf(file_ptr2, "M/N, Size, Fastest Algorithm, M, N, Mean Time(sec), Standard Deviation(sec), xFaster Combinatorial, xFaster Glynn, xFaster Ryser, Speed Combinatorial(sec), Speed Glynn(sec), Speed Ryser(sec) \n") < 0)
    {
        perror("Error occurred!");
        fclose(file_ptr2);
        return -1;
    }

    /* Generate the binomial coefficient table. */
    int64_t C[20][20];
    void bin_coeff();
    bin_coeff(20, 20, C);

    /* Store the maximum desired matrix dimensions from the input file in M and N. */

    int64_t M = m;
    int64_t N = n;
    double time_spent_on_comb[128];
    double time_spent_on_glynn[128];
    double time_spent_on_ryser[128];

    /* Specify the matrix dimension to solve. Matrices of sizes 2x2 up to MxN (specified by the user) must be generated. Remember m<=n so for each m we generate all possible n values before moving onto the next m. */

    for (int64_t m = 2; m <= M; m++)
    {
        for (int64_t n = m; n <= N; n++)
        {
            printf("Filling a %" PRId64 "-by-%" PRId64 " matrix with values and preparing to solve the permanent.\n", m, n);

            double matrix[256];
            double counter = 1.0;

            /* Populate the matrix of a given size with "dummy" values. These values will just go from 1 to mxn. */

            for (int64_t i = 0; i < m; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    matrix[i * n + j] = counter;
                    counter += 1.0;
                }
            }

            /* Solve the permanent using each algorithm the number of times specified by the user. */

            printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Combinatorial algorithm %" PRId64 " times.\n", m, n, r);

            double progress = 0.0;
            int64_t bar_width = 70;
            for (int64_t i = 0; i < r; i++)
            {
                clock_t begin_1 = clock(); // Time how long it takes to solve
                combinatorial((double *)matrix, m, n);
                clock_t end_1 = clock();
                time_spent_on_comb[i] = (double)(end_1 - begin_1) / CLOCKS_PER_SEC;

                /* Print progress bar. Adapted from C++ code found here: https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf */

                printf("[");
                progress += ((double)1 / (double)r);
                int pos = bar_width * progress;
                for (int64_t j = 0; j < bar_width; j++)
                {
                    if (j < pos)
                    printf("#");
                    else if (j == pos)
                    printf("#");
                    else
                    printf(" ");
                }
                printf("] %.2f, %%\r", (progress * (double)100));
            }
            printf("\n\n");

            printf("Combinatorial solution: %f \n\n", combinatorial((double *)matrix, m, n));
            
            printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Glynn algorithm %" PRId64 " times.\n", m, n, r);

            progress = 0.0;
            bar_width = 70;
            for (int64_t i = 0; i < r; i++)
            {
                clock_t begin_2 = clock();
                glynn((double *)matrix, m, n);
                clock_t end_2 = clock();
                time_spent_on_glynn[i] = (double)(end_2 - begin_2) / CLOCKS_PER_SEC;

                /* Print progress bar. */

                printf("[");
                progress += ((double)1 / (double)r);
                int pos = bar_width * progress;
                for (int64_t j = 0; j < bar_width; j++)
                {
                    if (j < pos)
                    printf("#");
                    else if (j == pos)
                    printf("#");
                    else
                    printf(" ");
                }
                printf("] %.2f, %%\r", (progress * (double)100));
            }
            printf("\n\n");

            printf("Glynn solution: %f \n\n", glynn((double *)matrix, m, n));
            
            printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Ryser algorithm %" PRId64 " times.\n", m, n, r);

            progress = 0.0;
            bar_width = 70;
            for (int64_t i = 0; i < r; i++)
            {
                clock_t begin_3 = clock();
                ryser((double *)matrix, m, n, (int64_t *)C);
                clock_t end_3 = clock();
                time_spent_on_ryser[i] = (double)(end_3 - begin_3) / CLOCKS_PER_SEC;

                /* Print progress bar. */

                printf("[");
                progress += ((double)1 / (double)r);
                int pos = bar_width * progress;
                for (int64_t j = 0; j < bar_width; j++)
                {
                    if (j < pos)
                    printf("#");
                    else if (j == pos)
                    printf("#");
                    else
                    printf(" ");
                }
                printf("] %.2f, %%\r", (progress * (double)100));
            }
            printf("\n\n");

            printf("Ryser solution: %f \n\n", ryser((double *)matrix, m, n, (int64_t *)C));

            double mean_time_comb = 0.0;
            double mean_time_glynn = 0.0;
            double mean_time_ryser = 0.0;

            /* Calculate the mean and standard deviation for the runtime of each algorithm. */

            for (int64_t i = 0; i < r; i++)
            {
                mean_time_comb += time_spent_on_comb[i]; // Sum up all of the time values
                mean_time_glynn += time_spent_on_glynn[i];
                mean_time_ryser += time_spent_on_ryser[i];
            }

            mean_time_comb = (double)mean_time_comb / (double)r; // Divide by the total number of runs done
            mean_time_glynn = (double)mean_time_glynn / (double)r;
            mean_time_ryser = (double)mean_time_ryser / (double)r;

            double sum_num_minus_mean_comb = 0.0;
            double sum_num_minus_mean_glynn = 0.0;
            double sum_num_minus_mean_ryser = 0.0;

            for (int i = 0; i < r; i++)
            {
                /* Sum all of the (values for runtime - mean runtime for each algorithm). */

                sum_num_minus_mean_comb += pow(time_spent_on_comb[i] - mean_time_comb, 2.0); 
                sum_num_minus_mean_glynn += pow(time_spent_on_glynn[i] - mean_time_glynn, 2.0);
                sum_num_minus_mean_ryser += pow(time_spent_on_ryser[i] - mean_time_ryser, 2.0);
            }

            /* Calculate the standard deviation for runtime of each algorithm. */

            double over_N = (double)1 / (double)r;
            double st_dev_comb = sqrt(over_N * sum_num_minus_mean_comb);
            double st_dev_glynn = sqrt(over_N * sum_num_minus_mean_glynn);
            double st_dev_ryser = sqrt(over_N * sum_num_minus_mean_ryser);
            
            /* Write all of the important information to the output file. */

            char s[] = "Fastest!";

            enum winning_algorithm alg;
            if (mean_time_comb <= mean_time_ryser && mean_time_comb <= mean_time_glynn)
            {
                alg = COMBINATORIAL;
            }
            else if (mean_time_glynn <= mean_time_ryser)
            {
                alg = GLYNN;
            }
            else 
            {
                alg = RYSER;
            }
            switch (alg)
            {
                case COMBINATORIAL:
                    if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Combinatorial,     %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %s, %.4fx, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_comb, st_dev_comb, s, (double)mean_time_glynn/(double)mean_time_comb, (double)mean_time_ryser/(double)mean_time_comb, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr2);
                        return -1;
                    }
                    break;
                case GLYNN:
                    if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Glynn,             %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %.4fx, %s, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_glynn, st_dev_glynn, (double)mean_time_comb/(double)mean_time_glynn, s, (double)mean_time_ryser/(double)mean_time_glynn, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr2);
                        return -1;
                    }
                    break;
                case RYSER:
                    if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Ryser,             %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %.4fx, %.4fx, %s, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_ryser, st_dev_ryser, (double)mean_time_comb/(double)mean_time_ryser, (double)mean_time_glynn/(double)mean_time_ryser, s, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr2);
                        return -1;
                    }
                    break;
            }
        }
    }
    printf("Fastest algorithm for solving the permanent of all matrices of sizes 2x2 up to %" PRIu64 "x%" PRIu64 " written successfully to file %s.\n", M, N, filename2);
    fclose(file_ptr2);
}
