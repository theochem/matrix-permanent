#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "permanent.h"


#define HEADER_FILE "permanent/tuning.h"

#define CSV_FILE    "fast_permanent.csv"

#define NUM_REPEATS 3

#define MAX_MATRIX  20

#define FASTEST     "Fastest!"

#define HEADER      "M/N, Size, Fastest Algorithm, M, N, Mean Time (sec), Standard Deviation (sec), Faster Combinatorial, Faster Glynn, Faster Ryser, Speed Combinatorial (sec), Speed Glynn (sec), Speed Ryser (sec)\n"

#define LINE_COMBN  "%.4e, %ld, Combn, %ld, %ld, %.10e, %.10e, %s, %.4ex, %.4ex, %.10e +- %.10e, %.10e +- %.10e, %.10e +- %.10e\n"

#define LINE_GLYNN  "%.4e, %ld, Glynn, %ld, %ld, %.10e, %.10e, %.4ex, %s, %.4ex, %.10e +- %.10e, %.10e +- %.10e, %.10e +- %.10e\n"

#define LINE_RYSER  "%.4e, %ld, Ryser, %ld, %ld, %.10e, %.10e, %.4ex, %.4ex, %s, %.10e +- %.10e, %.10e +- %.10e, %.10e +- %.10e\n"


double randl(double min, double max)
{
    srand((unsigned int) time(NULL));
    return ((double)rand() / (double)RAND_MAX) * (max - min) + min;
}


int main(void)
{
    /* Compute the tuning parameters */

    double param_1 = 3.14159;
    int64_t param_2 = 8192;

    /* Open a file for writing to */

    FILE *file_ptr = fopen(CSV_FILE, "w");
    if (file_ptr == NULL)
    {
        perror("Cannot open file!");
        return -1;
    }

    /* Print CSV headers */

    if (fprintf(file_ptr, HEADER) < 0)
    {
        perror("Error occurred!");
        fclose(file_ptr);
        return -1;
    }

    /* Time the efficiency of the algorithms for different size matrices. */

    int64_t i, m, n;
    double mean_time_combn, mean_time_glynn, mean_time_ryser;
    double sum_num_minus_mean_combn, sum_num_minus_mean_glynn, sum_num_minus_mean_ryser;
    double over_N, st_dev_combn, st_dev_glynn, st_dev_ryser;
    clock_t begin, end;
    double time_spent_on_combn[128];
    double time_spent_on_glynn[128];
    double time_spent_on_ryser[128];
    double randArray[4096];

    /* Random binary matrix for testing. */

    for (i = 0; i < 4096; i++)
    {
        randArray[i] = randl(-1.0, 1.0);
    }

    /* Iterate over number of rows and number of columns. */

    for (m = 2; m <= MAX_MATRIX; m++)
    {
        for (n = m; n <= MAX_MATRIX; n++)
        {
            /* Solve the permanent using each algorithm the number of times specified in NUM_REPEATS. */

            for (i = 0; i < NUM_REPEATS; i++)
            {
                begin = clock();
                combinatoric(m, n, randArray);
                end = clock();
                time_spent_on_combn[i] = (double)(end - begin) / (double)CLOCKS_PER_SEC;

                if (m == n)
                {
                    // begin = clock();
                    // glynn(m, n, randArray);
                    // end = clock();
                    // time_spent_on_glynn[i] = (double)(end - begin) / (double)CLOCKS_PER_SEC;
                    time_spent_on_glynn[i] = 1.0e9;

                    begin = clock();
                    ryser(m, n, randArray);
                    end = clock();
                    time_spent_on_ryser[i] = (double)(end - begin) / (double)CLOCKS_PER_SEC;
                    // time_spent_on_ryser[i] = 1.0e9;
                }
                else
                {
                    // begin = clock();
                    // glynn_rectangle(m, n, randArray);
                    // end = clock();
                    // time_spent_on_glynn[i] = (double)(end - begin) / (double)CLOCKS_PER_SEC;
                    time_spent_on_glynn[i] = 1.0e9;

                    begin = clock();
                    ryser_rectangle(m, n, randArray);
                    end = clock();
                    time_spent_on_ryser[i] = (double)(end - begin) / (double)CLOCKS_PER_SEC;
                    // time_spent_on_ryser[i] = 1.0e9;
                }
            }

            /* Calculate the mean and standard deviation for the runtime of each algorithm. */

            mean_time_combn = 0.0;
            mean_time_glynn = 0.0;
            mean_time_ryser = 0.0;

            for (i = 0; i < NUM_REPEATS; i++)
            {
                mean_time_combn += time_spent_on_combn[i];
                mean_time_glynn += time_spent_on_glynn[i];
                mean_time_ryser += time_spent_on_ryser[i];
            }

            mean_time_combn = (double)mean_time_combn / (double)NUM_REPEATS;
            mean_time_glynn = (double)mean_time_glynn / (double)NUM_REPEATS;
            mean_time_ryser = (double)mean_time_ryser / (double)NUM_REPEATS;

            sum_num_minus_mean_combn = 0.0;
            sum_num_minus_mean_glynn = 0.0;
            sum_num_minus_mean_ryser = 0.0;

            /* Sum all of the values for (runtime - mean runtime) for each algorithm. */

            for (i = 0; i < NUM_REPEATS; i++)
            {
                sum_num_minus_mean_combn += pow(time_spent_on_combn[i] - mean_time_combn, 2.0);
                sum_num_minus_mean_glynn += pow(time_spent_on_glynn[i] - mean_time_glynn, 2.0);
                sum_num_minus_mean_ryser += pow(time_spent_on_ryser[i] - mean_time_ryser, 2.0);
            }

            /* Calculate the standard deviation for runtime of each algorithm. */

            over_N = (double)1.0 / (double)NUM_REPEATS;
            st_dev_combn = sqrt(over_N * sum_num_minus_mean_combn);
            st_dev_glynn = sqrt(over_N * sum_num_minus_mean_glynn);
            st_dev_ryser = sqrt(over_N * sum_num_minus_mean_ryser);

            /* Write line depending on winning algorithm */

            if (mean_time_combn <= mean_time_ryser && mean_time_combn <= mean_time_glynn)
            {
                if (fprintf(file_ptr, LINE_COMBN, (double)m / (double)n, n, m, n, mean_time_combn, st_dev_combn, FASTEST, (double)mean_time_glynn / (double)mean_time_combn, (double)mean_time_ryser / (double)mean_time_combn, mean_time_combn, st_dev_combn, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                {
                    perror("Error occurred!");
                    fclose(file_ptr);
                    return -1;
                }
            }
            else if (mean_time_glynn <= mean_time_ryser)
            {
                if (fprintf(file_ptr, LINE_GLYNN, (double)m / (double)n, n, m, n, mean_time_glynn, st_dev_glynn, (double)mean_time_combn / (double)mean_time_glynn, FASTEST, (double)mean_time_ryser / (double)mean_time_glynn, mean_time_combn, st_dev_combn, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                {
                    perror("Error occurred!");
                    fclose(file_ptr);
                    return -1;
                }
            }
            else
            {
                if (fprintf(file_ptr, LINE_RYSER, (double)m / (double)n, n, m, n, mean_time_ryser, st_dev_ryser, (double)mean_time_combn / (double)mean_time_ryser, (double)mean_time_glynn / (double)mean_time_ryser, FASTEST, mean_time_combn, st_dev_combn, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
                {
                    perror("Error occurred!");
                    fclose(file_ptr);
                    return -1;
                }
            }
        }
    }

    /* Close written file */

    fclose(file_ptr);

    /* Write a header file with constants defined as macros */

    file_ptr = fopen(HEADER_FILE, "w");
    if (file_ptr == NULL)
    {
        perror("Cannot open file!");
        return -1;
    }
    if (fprintf(file_ptr, "#ifndef PERMANENT_TUNING_H\n#define PERMANENT_TUNING_H\n\n\n#define PARAM_1 %.9f\n#define PARAM_2 %ld\n\n\n#endif /* PERMANENT_TUNING_H */\n", param_1, param_2) < 0)
    {
        perror("Error occurred!");
        fclose(file_ptr);
        return -1;
    }

    /* Close written file */

    fclose(file_ptr);

    /* Exit successfully */

    return 0;
}
