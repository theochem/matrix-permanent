#include <cstdlib>
#include <cmath>
#include <ctime>
#include <random>
#include <algorithm>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "permanent.h"
#include "svm.h"


#ifdef RUN_TUNING

constexpr char CSV_FILE[] = "src/tuning.csv";

constexpr char CSV_HEADER[] = "M,N,Combn,Glynn,Ryser";

constexpr std::size_t NUM_TRIALS = 5;

constexpr std::size_t MAX_ROWS = 10;

constexpr std::size_t MAX_COLS = 10;

constexpr std::size_t DATA_POINTS = MAX_ROWS * MAX_COLS;



constexpr double TOLERANCE = 0.0001;

#else

constexpr char HEADER_FILE[] = "src/tuning.h";

constexpr double DEFAULT_PARAM_1 = 1.0;

constexpr double DEFAULT_PARAM_2 = 1.0;

constexpr double DEFAULT_PARAM_3 = 1.0;

constexpr double DEFAULT_PARAM_4 = 3.0;

#endif

/* Function to convert the first n rows of a 2D array into a vector of vectors */
std::vector<std::vector<double>> arrayToVector(int (*array)[2], int rows) {
    std::vector<std::vector<double>> result;
    for (int i = 0; i < rows; ++i) {
        std::vector<double> row;
        row.push_back((double)array[i][0]);
        row.push_back((double)array[i][1]);
        result.push_back(row);
    }

    return result;
}


/* Function to split the data into train and test sets */
void trainTestSplit(const std::vector<std::vector<double>>& data, 
                    std::vector<std::vector<double>>& trainData, 
                    std::vector<std::vector<double>>& testData,
                    double trainRatio) {
    /* Calculate the number of samples for the training set */
    int numTrainSamples = static_cast<int>(data.size() * trainRatio);

    /* Shuffle the data */
    std::vector<std::vector<double>> shuffledData = data;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(shuffledData.begin(), shuffledData.end(), gen);

    /* Split the shuffled data into train and test sets */
    trainData.assign(shuffledData.begin(), shuffledData.begin() + numTrainSamples);
    testData.assign(shuffledData.begin() + numTrainSamples, shuffledData.end());
}


#ifdef RUN_TUNING


int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    /* Declare array to hold tuning values */

    int metrics_combn[DATA_POINTS][2];
    int metrics_glynn[DATA_POINTS][2];
    int metrics_ryser[DATA_POINTS][2];

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
    double soln_ryser;

    std::clock_t begin, end;

    double time_combn[128];
    double time_glynn[128];
    double time_ryser[128];

    double mean_combn, mean_glynn, mean_ryser;

    int fastest;
    int counter_ryser =  0;
    int counter_combn =  0;
    int counter_glynn =  0;

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
                    time_combn[i] = (double)(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn(m, n, array);
                    end = std::clock();
                    time_glynn[i] = (double)(end - begin);

                    begin = std::clock();
                    soln_ryser = ryser(m, n, array);
                    end = std::clock();
                    time_ryser[i] = (double)(end - begin);

                    if (std::fabs(soln_combn - soln_glynn) / soln_combn > TOLERANCE) {
                        std::cerr << "Bad permanent values:"
                            << "\nCombn: " << soln_combn
                            << "\nGlynn: " << soln_glynn
                            << "\nRyser: " << soln_ryser << std::endl;
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
                    time_combn[i] = (double)(end - begin);

                    begin = std::clock();
                    soln_glynn = glynn_rectangular(m, n, array);
                    end = std::clock();
                    time_glynn[i] = (double)(end - begin);

                    time_ryser[i] = 1000000;

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

            /* Calculate the mean for the runtime of each algorithm. */

            mean_combn = 0.0;
            mean_glynn = 0.0;
            mean_ryser = 0.0;

            for (i = 0; i != NUM_TRIALS; ++i)
            {
                mean_combn += time_combn[i];
                mean_glynn += time_glynn[i];
                mean_ryser += time_ryser[i];
            }

            mean_combn = (double)mean_combn / (double)NUM_TRIALS;
            mean_glynn = (double)mean_glynn / (double)NUM_TRIALS;
            mean_ryser = (double)mean_ryser / (double)NUM_TRIALS;

            /* Find the fastest algorithm */

            if (mean_ryser <= mean_glynn) {
                metrics_ryser[counter_ryser][0] = m;
                metrics_ryser[counter_ryser][1] = n;
                counter_ryser += 1;
            } else if (mean_combn <= mean_glynn) {
                metrics_combn[counter_combn][0] = m;
                metrics_combn[counter_combn][1] = n;
                counter_combn += 1;
            } else {
                metrics_glynn[counter_glynn][0] = m;
                metrics_glynn[counter_glynn][1] = n;
                counter_glynn += 1;
            }

            /* Write line */

            std::cout << std::setw(3) << m << ',' << std::setw(3) << n << ','
                      << std::scientific << mean_combn << ',' << mean_glynn << ','
                      << std::scientific << mean_ryser << std::endl;

            csv_file << std::setw(3) << m << ',' << std::setw(3) << n << ','
                     << std::scientific << mean_combn << ',' << mean_glynn << ','
                     << std::scientific << mean_ryser << '\n';

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

    double DEFAULT_PARAM_4 = metrics_ryser[counter_ryser - 1][1];

    std::cout << std::setw(3) << DEFAULT_PARAM_4 << '\n';

    /* Prepare data for hard margin SVM boundary evaluation. */
    std::vector<std::vector<double>> clean_combn = arrayToVector(metrics_combn, counter_combn);
    std::vector<std::vector<double>> clean_glynn = arrayToVector(metrics_glynn, counter_glynn);

    /* Split into train/test for SVM */
    double trainRatio = 0.9; 
    std::vector<std::vector<double>> trainDataCombn;
    std::vector<std::vector<double>> testDataCombn;
    std::vector<std::vector<double>> trainDataGlynn;
    std::vector<std::vector<double>> testDataGlynn;

    /* Perform train-test split */
    trainTestSplit(clean_combn, trainDataCombn, testDataCombn, trainRatio);
    trainTestSplit(clean_glynn, trainDataGlynn, testDataGlynn, trainRatio);


    // /* Create train/test split */

    // std::vector<std::vector<double> > class1_data;
    // std::vector<std::vector<double>> class2_data;


    // std::vector<std::vector<double> > train_class1_data;
    // std::vector<std::vector<double>> train_class2_data;

    // (2.3) Training for SVM
    constexpr double lr = 0.01;
    size_t D = counter_combn + counter_glynn;
    HardMargin_SVM svm;
    svm.train(trainDataCombn, trainDataGlynn, D, lr);

    // // (3.1) Get Test Data for class 1
    // std::string test_class1_dir;
    // std::vector<std::string> test_class1_paths;
    // std::vector<std::vector<double>> test_class1_data;
    // /*****************************************************/
    // test_class1_dir = "datasets/" + vm["dataset"].as<std::string>() + "/" + vm["test_class1_dir"].as<std::string>();
    // test_class1_paths = Get_Paths(test_class1_dir);
    // test_class1_data = Get_Data(test_class1_paths, vm["nd"].as<size_t>());

    // // (3.2) Get Test Data for class 2
    // std::string test_class2_dir;
    // std::vector<std::string> test_class2_paths;
    // std::vector<std::vector<double>> test_class2_data;
    // /*****************************************************/
    // test_class2_dir = "datasets/" + vm["dataset"].as<std::string>() + "/" + vm["test_class2_dir"].as<std::string>();
    // test_class2_paths = Get_Paths(test_class2_dir);
    // test_class2_data = Get_Data(test_class2_paths, vm["nd"].as<size_t>());

    // (3.3) Test for SVM
    svm.test(testDataCombn, testDataGlynn);
    std::cout << "///////////////////////// Test /////////////////////////" << std::endl;
    std::cout << "accuracy: " << svm.accuracy << " (" << svm.correct_c1 + svm.correct_c2 << "/" << testDataCombn.size() + testDataGlynn.size() << ")" << std::endl;
    std::cout << "////////////////////////////////////////////////////////" << std::endl;




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
    header_file << "#define PARAM_4 " << std::scientific << DEFAULT_PARAM_4 << '\n';
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
