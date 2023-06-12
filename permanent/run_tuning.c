#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <tgmath.h>


#include "permanent.h"


enum winning_algorithm
{
    COMBINATORIAL,
    RYSER,
    GLYNN
};


int main()
{
    /* Compute the tuning parameters */
    double param_1 = 3.14159;
    int64_t param_2 = 8192;

    /* Open a file for writing to */

    const char filename[] = "fast_permanent.csv";
    FILE *file_ptr = NULL;

    file_ptr = fopen(filename, "w");
    if (file_ptr == NULL)
    {
        perror("Cannot open file!");
        return -1;
    }
    
    if (fprintf(file_ptr, "M/N, Size, Fastest Algorithm, M, N, Mean Time(sec), Standard Deviation(sec), xFaster Combinatorial, xFaster Glynn, xFaster Ryser, Speed Combinatorial(sec), Speed Glynn(sec), Speed Ryser(sec) \n") < 0)
    {
        perror("Error occurred!");
        fclose(file_ptr);
        return -1;
    }

    /* Time the efficiency of the algorithms for different size matrices. */
    double time_spent_on_comb[128];
    double time_spent_on_glynn[128];
    double time_spent_on_ryser[128];
    double time_c_2[128];
    double time_g_2[128];
    double time_r_2[128];
    int NUM_REPEATS = 3;
    int MAX_MATRIX = 20;

    for (int64_t m = 2; m <= MAX_MATRIX; m++)
    {
        for (int64_t n = m; n <= MAX_MATRIX; n++)
        {
            /* Populate the matrix of a given size randomly with 0's and 1's for testing. */
            double *randArray;
            double size = (double)(n * m);
            randArray = (double *)malloc(size * sizeof(double));
            if (randArray == NULL)
            {
                printf("Failed to allocate memory!\n");
                return -1;
            }

            /* Random binary matrix for testing. */
            for (int64_t i = 0; i < m; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    randArray[i * n + j ]= (double)(rand() % 2);
                }
            }

            /* Random 0, 2 matrix for testing. */
            double *randArray2;
            randArray2 = (double *)malloc(size * sizeof(double));
            if (randArray2 == NULL)
            {
                printf("Failed to allocate memory!\n");
                return -1;
            }
            for (int64_t i = 0; i < m; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    randArray2[i * n + j] = (double)(rand() % 2) * 2.0;
                }
            }

            /* Random -1, 1 matrix for testing. */
            double *randArray3;
            randArray3 = (double *)malloc(size * sizeof(double));
            if (randArray3 == NULL)
            {
                printf("Failed to allocate memory!\n");
                return -1;
            }
            for (int64_t i = 0; i < m; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    randArray3[i * n + j] = -1.0 + (double)(rand() % 2) * 2.0;
                }
            }

            /* Random -2, 2 matrix for testing. */
            double *randArray4;
            randArray4 = (double *)malloc(size * sizeof(double));
            if (randArray4 == NULL)
            {
                printf("Failed to allocate memory!\n");
                return -1;
            }
            for (int64_t i = 0; i < m; i++)
            {
                for (int64_t j = 0; j < n; j++)
                {
                    randArray4[i * n + j] = -2.0 + (double)(rand() % 2) * 4.0;
                }
            }

            /* Solve the permanent using each algorithm the number of times specified in NUM_REPEATS. */
            for (int64_t i = 0; i < NUM_REPEATS; i++)
            {
                //double squareness = (double)m / (double)n ;
                double current_size = (double)m ;
                double max_size = 6.0;
                double comparison_value = 0.6;
                if ((current_size <= max_size) || (current_size < comparison_value))
                {
                    clock_t begin_1 = clock(); // Time how long it takes to solve
                    combinatoric(m, n, randArray);
                    clock_t end_1 = clock();
                    time_spent_on_comb[i] = (double)(end_1 - begin_1) / CLOCKS_PER_SEC;

                    clock_t begin_12 = clock(); // Time how long it takes to solve
                    combinatoric(m, n, randArray2);
                    clock_t end_12 = clock();
                    time_c_2[i] = (double)(end_12 - begin_12) / CLOCKS_PER_SEC;
                }
                else
                {
                    time_spent_on_comb[i] = 100.0;
                    time_c_2[i] = 100.0;
                }

                if (m == n)
                {
                    clock_t begin_2 = clock();
                    glynn(m, n, randArray);
                    clock_t end_2 = clock();
                    time_spent_on_glynn[i] = (double)(end_2 - begin_2) / CLOCKS_PER_SEC;

                    clock_t begin_22 = clock(); // Time how long it takes to solve
                    glynn(m, n, randArray2);
                    clock_t end_22 = clock();
                    time_g_2[i] = (double)(end_22 - begin_22) / CLOCKS_PER_SEC;

                    clock_t begin_3 = clock();
                    ryser(m, n, randArray);
                    clock_t end_3 = clock();
                    time_spent_on_ryser[i] = (double)(end_3 - begin_3) / CLOCKS_PER_SEC;

                    clock_t begin_32 = clock(); // Time how long it takes to solve
                    ryser(m, n, randArray2);
                    clock_t end_32 = clock();
                    time_r_2[i] = (double)(end_32 - begin_32) / CLOCKS_PER_SEC;
                }

                else if (m != n)
                {
                    clock_t begin_2 = clock();
                    glynn_rectangle(m, n, randArray);
                    clock_t end_2 = clock();
                    time_spent_on_glynn[i] = (double)(end_2 - begin_2) / CLOCKS_PER_SEC;

                    clock_t begin_22 = clock(); // Time how long it takes to solve
                    glynn_rectangle(m, n, randArray2);
                    clock_t end_22 = clock();
                    time_g_2[i] = (double)(end_22 - begin_22) / CLOCKS_PER_SEC;

                    if (current_size <= max_size)
                    {
                        clock_t begin_3 = clock();
                        ryser_rectangle(m, n, randArray);
                        clock_t end_3 = clock();
                        time_spent_on_ryser[i] = (double)(end_3 - begin_3) / CLOCKS_PER_SEC;

                        clock_t begin_32 = clock(); // Time how long it takes to solve
                        ryser_rectangle(m, n, randArray2);
                        clock_t end_32 = clock();
                        time_r_2[i] = (double)(end_32 - begin_32) / CLOCKS_PER_SEC;
                    }
                    else
                    {
                        time_spent_on_ryser[i] = 100.0;
                        time_r_2[i] = 100.0;
                    }
                }
            }

            double mean_time_comb = 0.0;
            double mean_time_glynn = 0.0;
            double mean_time_ryser = 0.0;

            /* Calculate the mean and standard deviation for the runtime of each algorithm. */
            for (int64_t i = 0; i < NUM_REPEATS; i++)
            {
                mean_time_comb += time_spent_on_comb[i]; // Sum up all of the time values
                mean_time_glynn += time_spent_on_glynn[i];
                mean_time_ryser += time_spent_on_ryser[i];

                mean_time_comb += time_c_2[i];
                mean_time_glynn += time_g_2[i];
                mean_time_ryser += time_r_2[i];
            }
            
            mean_time_comb = (double)mean_time_comb / (double)NUM_REPEATS; // Divide by the total number of runs done
            mean_time_glynn = (double)mean_time_glynn / (double)NUM_REPEATS;
            mean_time_ryser = (double)mean_time_ryser / (double)NUM_REPEATS;

            double sum_num_minus_mean_comb = 0.0;
            double sum_num_minus_mean_glynn = 0.0;
            double sum_num_minus_mean_ryser = 0.0;

            for (int i = 0; i < NUM_REPEATS; i++)
            {
                /* Sum all of the (values for runtime - mean runtime for each algorithm). */
                sum_num_minus_mean_comb += pow(time_spent_on_comb[i] - mean_time_comb, 2.0); 
                sum_num_minus_mean_glynn += pow(time_spent_on_glynn[i] - mean_time_glynn, 2.0);
                sum_num_minus_mean_ryser += pow(time_spent_on_ryser[i] - mean_time_ryser, 2.0);
            }

            /* Calculate the standard deviation for runtime of each algorithm. */
            double over_N = (double)1 / (double)NUM_REPEATS;
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
                    if (fprintf(file_ptr, "%.4f, %lld, Combinatorial,     %lld, %lld, %.10f, %.10f, %s, %.4fx, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_comb, st_dev_comb, s, (double)mean_time_glynn/(double)mean_time_comb, (double)mean_time_ryser/(double)mean_time_comb, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr);
                        return -1;
                    }
                    break;
                case GLYNN:
                    if (fprintf(file_ptr, "%.4f, %lld, Glynn,             %lld, %lld, %.10f, %.10f, %.4fx, %s, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_glynn, st_dev_glynn, (double)mean_time_comb/(double)mean_time_glynn, s, (double)mean_time_ryser/(double)mean_time_glynn, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr);
                        return -1;
                    }
                    break;
                case RYSER:
                    if (fprintf(file_ptr, "%.4f, %lld, Ryser,             %lld, %lld, %.10f, %.10f, %.4fx, %.4fx, %s, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_ryser, st_dev_ryser, (double)mean_time_comb/(double)mean_time_ryser, (double)mean_time_glynn/(double)mean_time_ryser, s, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                    {
                        perror("Error occurred!");
                        fclose(file_ptr);
                        return -1;
                    }
                    break;
            }
            //fprintf(file_ptr, "Successfully allocated memory for an array of %f elements of %zu bytes each at address %p using malloc.\n", size, sizeof(*randArray), (void *)randArray);
            //fprintf(file_ptr, "First element (uninitialized): %f.\n", randArray[0]);
            //fprintf(file_ptr, "Check address of void pointer: %p and double pointer: %p. \n", (void *)randArray, (double *)randArray);
            // if (m == n)
            // {
            //     fprintf(file_ptr, "Solution to Glynn: %f, Solution to Ryser: %f \n", glynn(m, n, randArray), ryser(m, n, randArray));
            // }
            // else if (m != n)
            // {
            //     fprintf(file_ptr, "Solution to Glynn: %f, Solution to Ryser: %f \n", glynn_rectangle(m, n, randArray), ryser_rectangle(m, n, randArray));
            // }
            
            free(randArray);
            free(randArray2);
            free(randArray3);
            free(randArray4);
        }
    }
    fclose(file_ptr);

    /* Print a header file to stdout with constants defined as macros */
    printf("#ifndef PERMANENT_TUNING_H\n");
    printf("#define PERMANENT_TUNING_H\n");
    printf("\n\n");
    printf("#define PARAM_1 %.9f\n", param_1);
    printf("#define PARAM_2 %lld\n", param_2);
    printf("\n\n");
    printf("#endif /* PERMANENT_TUNING_H */\n");

    /* Exit successfully */
    return 0;
}
