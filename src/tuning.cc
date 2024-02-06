#include <cstdlib>
#include <cmath>
#include <ctime>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "permanent.h"


#ifdef RUN_TUNING

constexpr char CSV_FILE[] = "src/tuning.csv";

constexpr char CSV_HEADER[] = "M,N,Combn,Glynn,Fastest";

constexpr std::size_t NUM_TRIALS = 3;

constexpr std::size_t MAX_ROWS = 10;

constexpr std::size_t MAX_COLS = 16;

constexpr double TOLERANCE = 0.0001;

#else

constexpr char HEADER_FILE[] = "src/tuning.h";

constexpr double DEFAULT_PARAM_1 = 1.0;

constexpr double DEFAULT_PARAM_2 = 1.0;

constexpr double DEFAULT_PARAM_3 = 1.0;

#endif


#ifdef RUN_TUNING


int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    /* Open CSV file */

    std::ofstream csv_file;

    csv_file.open(CSV_FILE);

    if (csv_file.fail())
    {
        std::cerr << "Cannot open CSV file" << std::endl;
        return 2;
    }

    /* Set format options for printing */

    std::cout.precision(9);
    std::cerr.precision(9);
    csv_file.precision(9);

    std::cout.fill(' ');
    std::cerr.fill(' ');
    csv_file.fill(' ');

    /* Print CSV headers */

    std::cout << CSV_HEADER << std::endl;

    csv_file << CSV_HEADER << '\n';

    if (csv_file.fail())
    {
        std::cerr << "Error writing to CSV file" << std::endl;
        csv_file.close();
        return 2;
    }

    /* Time efficiency of algorithms for different size matrices */

    std::size_t i, m, n;

    double soln_combn;
    double soln_glynn;

    std::clock_t begin, end;

    double time_combn;
    double time_glynn;

    int fastest;

    double array[MAX_ROWS * MAX_COLS];

    /* Random binary matrix for testing. */

    std::srand(1234789789U);
    for (i = 0; i < MAX_ROWS * MAX_COLS; ++i)
    {
        array[i] = (std::rand() / (double)RAND_MAX - 0.5) * 2;
    }

    /* Iterate over number of rows and number of columns. */

    for (m = 2; m <= MAX_ROWS; ++m)
    {
        for (n = m; n <= MAX_COLS; ++n)
        {
            /* Solve the permanent using each algorithm NUM_TRIALS number of times. */

            if (m == n)
            {
                for (i = 0; i != NUM_TRIALS; ++i)
                {

                    begin = std::clock();
                    soln_combn = combinatoric(m, n, array);
                    end = std::clock();
                    time_combn = (double)(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn(m, n, array);
                    end = std::clock();
                    time_glynn = (double)(end - begin);

                    if (std::fabs(soln_combn - soln_glynn) / soln_combn > TOLERANCE) {
                        std::cerr << "Bad permanent values:"
                            << "\nCombn: " << soln_combn
                            << "\nGlynn: " << soln_glynn << std::endl;
                        csv_file.close();
                        return 1;
                    }
                }
            }
            else
            {
                for (i = 0; i != NUM_TRIALS; ++i)
                {
                    begin = std::clock();
                    soln_combn = combinatoric_rectangular(m, n, array);
                    end = std::clock();
                    time_combn = (double)(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn_rectangular(m, n, array);
                    end = std::clock();
                    time_glynn = (double)(end - begin);

                    if (std::fabs(soln_combn - soln_glynn) / soln_combn > TOLERANCE) {
                        std::cerr << std::scientific
                            << "Bad permanent values:"
                            << "\nCombn: " << soln_combn
                            << "\nGlynn: " << soln_glynn << std::endl;
                        csv_file.close();
                        return 1;
                    }
                }
            }

            /* Find the fastest algorithm */

            if (time_combn <= time_glynn) {
                fastest = 0;
            } else {
                fastest = 1;
            }

            /* Write line */

            std::cout << std::setw(3) << m << ',' << std::setw(3) << n << ','
                      << std::scientific << time_combn << ',' << time_glynn << ','
                      << fastest << std::endl;

            csv_file << std::setw(3) << m << ',' << std::setw(3) << n << ','
                     << std::scientific << time_combn << ',' << time_glynn << ','
                     << fastest << '\n';

            if (csv_file.fail())
            {
                std::cerr << "Error writing to CSV file" << std::endl;
                csv_file.close();
                return 2;
            }
        }
    }

    /* Close CSV file */

    csv_file.close();

    /* Exit successfully */

    return 0;
}


#else


int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    /* Open header file */

    std::ofstream header_file;

    header_file.open(HEADER_FILE);

    if (header_file.fail())
    {
        std::cerr << "Cannot open header file" << std::endl;
        return 2;
    }

    /* Set format options for printing */

    header_file.precision(9);

    /* Write default header file */

    header_file << "#ifndef PERMANENT_TUNING_H\n";
    header_file << "#define PERMANENT_TUNING_H\n";
    header_file << "\n\n";
    header_file << "#define PARAM_1 " << std::scientific << DEFAULT_PARAM_1 << '\n';
    header_file << "#define PARAM_2 " << std::scientific << DEFAULT_PARAM_2 << '\n';
    header_file << "#define PARAM_3 " << std::scientific << DEFAULT_PARAM_3 << '\n';
    header_file << "\n\n";
    header_file << "#endif /* PERMANENT_TUNING_H */\n";

    if (header_file.fail()) {
        std::cerr << "Error writing to header file" << std::endl;
        header_file.close();
        return 2;
    }

    /* Close header file */

    header_file.close();

    /* Exit successfully */

    return 0;
}


#endif
