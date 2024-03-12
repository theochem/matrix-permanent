/* Copyright 2034 QC-Devs (GPLv3) */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "permanent.h"

#ifdef WITH_TUNING_FILE

#include "tuning.h"

#endif  // WITH_TUNING_FILE

#ifdef RUN_TUNING

constexpr char CSV_FILE[] = "src/tuning.csv";

constexpr char CSV_HEADER[] = "M/N,N,Combn,Glynn,Ryser,Fastest";

constexpr std::size_t NUM_TRIALS = 5;

constexpr std::size_t MAX_ROWS = 16;

constexpr std::size_t MAX_COLS = 16;

constexpr std::size_t DATA_POINTS = MAX_ROWS * MAX_COLS;

constexpr double TOLERANCE = 0.0001;

#else

constexpr char HEADER_FILE[] = "src/tuning.h";

constexpr double DEFAULT_PARAM_1 = -0.572098;

constexpr double DEFAULT_PARAM_2 = -22.014212;

constexpr double DEFAULT_PARAM_3 = 15.297940;

constexpr double DEFAULT_PARAM_4 = 3.0;

#endif  // RUN_TUNING

#ifdef RUN_TUNING

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    /* Open CSV file */

    std::ofstream csv_file;

    csv_file.open(CSV_FILE);

    if (csv_file.fail()) {
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

    if (csv_file.fail()) {
        std::cerr << "Error writing to CSV file" << std::endl;
        csv_file.close();
        return 2;
    }

    /* Time efficiency of algorithms for different size matrices */

    std::size_t i, m, n;

    std::clock_t begin, end;

    double soln_combn, soln_glynn, soln_ryser;
    double mean_combn, mean_glynn, mean_ryser;

    double time_combn[128];
    double time_glynn[128];
    double time_ryser[128];

    double array[MAX_ROWS * MAX_COLS];

    int fastest;

    /* Random binary matrix for testing. */

    std::srand(1234789789U);
    for (i = 0; i < MAX_ROWS * MAX_COLS; ++i) {
        array[i] = (std::rand() / static_cast<double>(RAND_MAX) - 0.5) * 2;
    }

    /* Iterate over number of rows and number of columns. */

    for (m = 2; m <= MAX_ROWS; ++m) {
        for (n = m; n <= MAX_COLS; ++n) {
            /* Solve the permanent using each algorithm NUM_TRIALS number of
             * times. */

            if (m == n) {
                for (i = 0; i != NUM_TRIALS; ++i) {
                    begin = std::clock();
                    soln_combn = combinatoric<double>(m, n, array);
                    end = std::clock();
                    time_combn[i] = static_cast<double>(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn<double>(m, n, array);
                    end = std::clock();
                    time_glynn[i] = static_cast<double>(end - begin);

                    begin = std::clock();
                    soln_ryser = ryser<double>(m, n, array);
                    end = std::clock();
                    time_ryser[i] = static_cast<double>(end - begin);

                    if ((std::fabs(soln_combn - soln_glynn) / soln_combn >
                         TOLERANCE) ||
                        (std::fabs(soln_combn - soln_ryser) / soln_ryser >
                         TOLERANCE)) {
                        std::cerr << "Bad permanent values:"
                                  << "\nCombn: " << soln_combn
                                  << "\nGlynn: " << soln_glynn
                                  << "\nRyser: " << soln_ryser << std::endl;
                        csv_file.close();
                        return 1;
                    }
                }
            } else {
                for (i = 0; i != NUM_TRIALS; ++i) {
                    begin = std::clock();
                    soln_combn = combinatoric_rectangular<double>(m, n, array);
                    end = std::clock();
                    time_combn[i] = static_cast<double>(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn_rectangular<double>(m, n, array);
                    end = std::clock();
                    time_glynn[i] = static_cast<double>(end - begin);

                    begin = std::clock();
                    soln_ryser = ryser_rectangular<double>(m, n, array);
                    end = std::clock();
                    time_ryser[i] = static_cast<double>(end - begin);

                    if ((std::fabs(soln_combn - soln_glynn) / soln_combn >
                         TOLERANCE) ||
                        (std::fabs(soln_combn - soln_ryser) / soln_ryser >
                         TOLERANCE)) {
                        std::cerr << std::scientific << "Bad permanent values:"
                                  << "\nCombn: " << soln_combn
                                  << "\nGlynn: " << soln_glynn
                                  << "\nRyser: " << soln_ryser << std::endl;
                        csv_file.close();
                        return 1;
                    }
                }
            }

            /* Calculate the mean for the runtime of each algorithm. */

            mean_combn = 0.0;
            mean_glynn = 0.0;
            mean_ryser = 0.0;

            for (i = 0; i != NUM_TRIALS; ++i) {
                mean_combn += time_combn[i];
                mean_glynn += time_glynn[i];
                mean_ryser += time_ryser[i];
            }

            mean_combn = mean_combn / NUM_TRIALS;
            mean_glynn = mean_glynn / NUM_TRIALS;
            mean_ryser = mean_ryser / NUM_TRIALS;

            /* Find the fastest algorithm */

            if (mean_ryser <= mean_glynn) {
                fastest = 0;
            } else if (mean_combn <= mean_glynn) {
                fastest = 1;
            } else {
                fastest = 2;
            }

            /* Write line */

            std::cout << std::scientific << static_cast<double>(m) / n << ','
                      << std::setw(2) << n << ',' << std::scientific
                      << mean_combn << ',' << mean_glynn << ','
                      << std::scientific << mean_ryser << ',' << fastest
                      << std::endl;

            csv_file << std::scientific << static_cast<double>(m) / n << ','
                     << std::setw(2) << n << ',' << std::scientific
                     << mean_combn << ',' << mean_glynn << ','
                     << std::scientific << mean_ryser << ',' << fastest << '\n';

            if (csv_file.fail()) {
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

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    /* Open header file */

    std::ofstream header_file;

    header_file.open(HEADER_FILE);

    if (header_file.fail()) {
        std::cerr << "Cannot open header file" << std::endl;
        return 2;
    }

    /* Set format options for printing */

    header_file.precision(9);

    /* Write default header file */

    header_file << "#ifndef TUNING_H_\n";
    header_file << "#define TUNING_H_\n";
    header_file << "\n\n";
    header_file << "constexpr double PARAM_1 = " << std::scientific
                << DEFAULT_PARAM_1 << ";\n";
    header_file << "constexpr double PARAM_2 = " << std::scientific
                << DEFAULT_PARAM_2 << ";\n";
    header_file << "constexpr double PARAM_3 = " << std::scientific
                << DEFAULT_PARAM_3 << ";\n";
    header_file << "constexpr double PARAM_4 = " << std::scientific
                << DEFAULT_PARAM_4 << ";\n";
    header_file << "\n\n";
    header_file << "#endif  // TUNING_H_\n";

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

#endif  // RUN_TUNING
