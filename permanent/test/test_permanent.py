import numpy as np
import permanent

def bin_coeff():
    C = np.empty(shape=(100, 100), dtype='object')
    for k in range(1, 101):
        C[0][k] = 0
    
    for n in range(0, 101):
        C[n][0] = 1

    for n in range(1, 101):
        for k in range(1, 101):
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k]
    
    return C


def test_2by2_comb():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_comb():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.combinatoric(matrix)==450

def test_4by4_comb():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.combinatoric(matrix)==55456

def test_2by3_comb():
    matrix = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.double)
    assert permanent.combinatoric(matrix)==58.0

def test_2by2_glynn():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.combinatoric(matrix)==10

def test_3by3_glynn():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix)==450

def test_4by4_glynn():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix)==55456

def test_2by2_ryser():
    matrix = np.arange(1, 5, dtype=np.double).reshape(2,2)
    assert permanent.ryser(matrix, bin_coeff())==10

def test_3by3_ryser():
    matrix = np.arange(1, 10, dtype=np.double).reshape(3,3)
    assert permanent.ryser(matrix, bin_coeff())==450

def test_4by4_ryser():
    matrix = np.arange(1, 17, dtype=np.double).reshape(4,4)
    assert permanent.ryser(matrix, bin_coeff())==55456

def test_2by3_ryser():
    matrix = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.double)
    assert permanent.ryser(matrix, bin_coeff())==58.0



#     ## Run a speed test!!

#     enum winning_algorithm
# {
#     COMBINATORIAL,
#     RYSER,
#     GLYNN
# };

#     static PyObject *permanent(PyObject *module, PyObject *object)
# {
#     /* Cast PyObject* to PyArrayObject*. */
#     PyArrayObject *matrix = (PyArrayObject *)PyArray_FromAny(object, NULL, 2, 2, NPY_ARRAY_ALIGNED, NULL);
#     double *ptr = (double *)PyArray_GETPTR2(matrix, 0, 0);

#     /* Open a file for writing to. */

#     const char filename2[] = "fast_permanent.csv";
#     FILE *file_ptr2 = NULL;

#     file_ptr2 = fopen(filename2, "w");
#     if (file_ptr2 == NULL)
#     {
#         perror("Oh no! Cannot open file!");
#         return -1;
#     }
    
#     if (fprintf(file_ptr2, "M/N, Size, Fastest Algorithm, M, N, Mean Time(sec), Standard Deviation(sec), xFaster Combinatorial, xFaster Glynn, xFaster Ryser, Speed Combinatorial(sec), Speed Glynn(sec), Speed Ryser(sec) \n") < 0)
#     {
#         perror("Error occurred!");
#         fclose(file_ptr2);
#         return -1;
#     }

#     // /* Generate the binomial coefficient table. */
#     // int64_t C[100][100];
#     // void bin_coeff();
#     // bin_coeff(100, 100, C);

#     /* Store the maximum desired matrix dimensions from the input file in M and N. */

#     int64_t m_rows = PyArray_DIMS(matrix)[0];
#     int64_t n_cols = PyArray_DIMS(matrix)[1];

#     double time_spent_on_comb[128];
#     double time_spent_on_glynn[128];
#     double time_spent_on_ryser[128];

#     /* Specify the matrix dimension to solve. Matrices of sizes 2x2 up to MxN (specified by the user) must be generated. Remember m<=n so for each m we generate all possible n values before moving onto the next m. */

#     for (int64_t m = 2; m <= m_rows; m++)
#     {
#         for (int64_t n = m; n <= n_cols; n++)
#         {
#             printf("Filling a %" PRId64 "-by-%" PRId64 " matrix with values and preparing to solve the permanent.\n", m, n);

#             double matrix[256];
#             double counter = 1.0;

#             /* Populate the matrix of a given size with "dummy" values. These values will just go from 1 to mxn. */

#             for (int64_t i = 0; i < m; i++)
#             {
#                 for (int64_t j = 0; j < n; j++)
#                 {
#                     matrix[i * n + j] = counter;
#                     counter += 1.0;
#                 }
#             }

#             /* Solve the permanent using each algorithm the number of times specified by the user. */

#             printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Combinatorial algorithm 10 times.\n", m, n);

#             double progress = 0.0;
#             int64_t bar_width = 70;
#             for (int64_t i = 0; i < 10; i++)
#             {
#                 clock_t begin_1 = clock(); // Time how long it takes to solve
#                 combinatorial(matrix);
#                 clock_t end_1 = clock();
#                 time_spent_on_comb[i] = (double)(end_1 - begin_1) / CLOCKS_PER_SEC;

#                 /* Print progress bar. Adapted from C++ code found here: https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf */

#                 printf("[");
#                 progress += ((double)1 / (double)10);
#                 int pos = bar_width * progress;
#                 for (int64_t j = 0; j < bar_width; j++)
#                 {
#                     if (j < pos)
#                     printf("#");
#                     else if (j == pos)
#                     printf("#");
#                     else
#                     printf(" ");
#                 }
#                 printf("] %.2f, %%\10", (progress * (double)100));
#             }
#             printf("\n\n");

#             printf("Combinatorial solution: %f \n\n", combinatorial((double *)matrix));
            
#             printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Glynn algorithm 10 times.\n", m, n);

#             progress = 0.0;
#             bar_width = 70;
#             for (int64_t i = 0; i < 10; i++)
#             {
#                 clock_t begin_2 = clock();
#                 glynn(matrix);
#                 clock_t end_2 = clock();
#                 time_spent_on_glynn[i] = (double)(end_2 - begin_2) / CLOCKS_PER_SEC;

#                 /* Print progress bar. */

#                 printf("[");
#                 progress += ((double)1 / (double)10);
#                 int pos = bar_width * progress;
#                 for (int64_t j = 0; j < bar_width; j++)
#                 {
#                     if (j < pos)
#                     printf("#");
#                     else if (j == pos)
#                     printf("#");
#                     else
#                     printf(" ");
#                 }
#                 printf("] %.2f, %%\10", (progress * (double)100));
#             }
#             printf("\n\n");

#             printf("Glynn solution: %f \n\n", glynn((double *)matrix);
            
#             printf("Solving the permanent of a %" PRId64 "-by-%" PRId64 " matrix using the Ryser algorithm 10 times.\n", m, n);

#             progress = 0.0;
#             bar_width = 70;
#             for (int64_t i = 0; i < 10; i++)
#             {
#                 clock_t begin_3 = clock();
#                 ryser(matrix, C);
#                 clock_t end_3 = clock();
#                 time_spent_on_ryser[i] = (double)(end_3 - begin_3) / CLOCKS_PER_SEC;

#                 /* Print progress bar. */

#                 printf("[");
#                 progress += ((double)1 / (double)10);
#                 int pos = bar_width * progress;
#                 for (int64_t j = 0; j < bar_width; j++)
#                 {
#                     if (j < pos)
#                     printf("#");
#                     else if (j == pos)
#                     printf("#");
#                     else
#                     printf(" ");
#                 }
#                 printf("] %.2f, %%\10", (progress * (double)100));
#             }
#             printf("\n\n");

#             printf("Ryser solution: %f \n\n", ryser((double *)matrix, (double *)C));

#             double mean_time_comb = 0.0;
#             double mean_time_glynn = 0.0;
#             double mean_time_ryser = 0.0;

#             /* Calculate the mean and standard deviation for the runtime of each algorithm. */

#             for (int64_t i = 0; i < 10; i++)
#             {
#                 mean_time_comb += time_spent_on_comb[i]; // Sum up all of the time values
#                 mean_time_glynn += time_spent_on_glynn[i];
#                 mean_time_ryser += time_spent_on_ryser[i];
#             }

#             mean_time_comb = (double)mean_time_comb / (double)10; // Divide by the total number of runs done
#             mean_time_glynn = (double)mean_time_glynn / (double)10;
#             mean_time_ryser = (double)mean_time_ryser / (double)10;

#             double sum_num_minus_mean_comb = 0.0;
#             double sum_num_minus_mean_glynn = 0.0;
#             double sum_num_minus_mean_ryser = 0.0;

#             for (int i = 0; i < 10; i++)
#             {
#                 /* Sum all of the (values for runtime - mean runtime for each algorithm). */

#                 sum_num_minus_mean_comb += pow(time_spent_on_comb[i] - mean_time_comb, 2.0); 
#                 sum_num_minus_mean_glynn += pow(time_spent_on_glynn[i] - mean_time_glynn, 2.0);
#                 sum_num_minus_mean_ryser += pow(time_spent_on_ryser[i] - mean_time_ryser, 2.0);
#             }

#             /* Calculate the standard deviation for runtime of each algorithm. */

#             double over_N = (double)1 / (double)10;
#             double st_dev_comb = sqrt(over_N * sum_num_minus_mean_comb);
#             double st_dev_glynn = sqrt(over_N * sum_num_minus_mean_glynn);
#             double st_dev_ryser = sqrt(over_N * sum_num_minus_mean_ryser);
            
#             /* Write all of the important information to the output file. */

#             char s[] = "Fastest!";

#             enum winning_algorithm alg;
#             if (mean_time_comb <= mean_time_ryser && mean_time_comb <= mean_time_glynn)
#             {
#                 alg = COMBINATORIAL;
#             }
#             else if (mean_time_glynn <= mean_time_ryser)
#             {
#                 alg = GLYNN;
#             }
#             else 
#             {
#                 alg = RYSER;
#             }
#             switch (alg)
#             {
#                 case COMBINATORIAL:
#                     if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Combinatorial,     %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %s, %.4fx, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_comb, st_dev_comb, s, (double)mean_time_glynn/(double)mean_time_comb, (double)mean_time_ryser/(double)mean_time_comb, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
#                     {
#                         perror("Error occurred!");
#                         fclose(file_ptr2);
#                         return -1;
#                     }
#                     break;
#                 case GLYNN:
#                     if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Glynn,             %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %.4fx, %s, %.4fx, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_glynn, st_dev_glynn, (double)mean_time_comb/(double)mean_time_glynn, s, (double)mean_time_ryser/(double)mean_time_glynn, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
#                     {
#                         perror("Error occurred!");
#                         fclose(file_ptr2);
#                         return -1;
#                     }
#                     break;
#                 case RYSER:
#                     if (fprintf(file_ptr2, "%.4f, %"PRIu64 ", Ryser,             %"PRIu64 ", %"PRIu64 ", %.10f, %.10f, %.4fx, %.4fx, %s, %.10f +- %.10f, %.10f +- %.10f, %.10f +- %.10f \n", (double)m/(double)n, n, m, n, mean_time_ryser, st_dev_ryser, (double)mean_time_comb/(double)mean_time_ryser, (double)mean_time_glynn/(double)mean_time_ryser, s, mean_time_comb, st_dev_comb, mean_time_glynn, st_dev_glynn, mean_time_ryser, st_dev_ryser) < 0)
#                     {
#                         perror("Error occurred!");
#                         fclose(file_ptr2);
#                         return -1;
#                     }
#                     break;
#             }
#         }
#     }
#     printf("Fastest algorithm for solving the permanent of all matrices of sizes 2x2 up to %" PRIu64 "x%" PRIu64 " written successfully to file %s.\n", m_rows, n_cols, filename2);
#     fclose(file_ptr2);
# };