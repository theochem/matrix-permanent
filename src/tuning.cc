/* Copyright 2024-2025 QC-Devs (GPLv3) */

#if defined(__GNUC__) || defined(__clang__)
#define _PERMANENT_ALWAYS_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER) && !defined(__clang__)
#define _PERMANENT_ALWAYS_INLINE __forceinline
#define __func__ __FUNCTION__
#else
#define _PERMANENT_ALWAYS_INLINE
#endif

#define _PERMANENT_DEFAULT_TUNING

#include <permanent.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace permanent {

constexpr char CSV_FILE[] = "src/tuning.csv";

constexpr char CSV_HEADER[] = "M/N,N,Combn,Glynn,Ryser,Fastest";

constexpr size_t NUM_TRIALS = 5;

constexpr size_t MAX_ROWS = 26;

constexpr size_t MAX_COLS = 26;

constexpr size_t DATA_POINTS = MAX_ROWS * MAX_COLS;

constexpr double TOLERANCE = 0.0001;

constexpr size_t RUN_COMBINATORIAL_UNTIL = 13;

namespace {

template <typename T>
inline _PERMANENT_ALWAYS_INLINE void ensure(T &value)
{
  asm volatile("" : "+r,m"(value) : : "memory");
}

template <typename T, typename I = void>
int generate_tuning_data(std::ofstream &csv_file, const size_t m, const size_t n, const T *array,
                         double *time_combn, double *time_glynn, double *time_ryser)
{
  // Solve the permanent using each algorithm NUM_TRIALS number of times
  result_t<T, I> soln_combn, soln_glynn, soln_ryser;
  std::clock_t begin, end;
  bool run_combinatorial = (m <= RUN_COMBINATORIAL_UNTIL);
  if (m == n) {
    for (size_t i = 0; i != NUM_TRIALS; ++i) {
      if (run_combinatorial) {
        begin = std::clock();
        soln_combn = permanent::combinatoric_square<T, I>(m, n, array);
        ensure(soln_combn);
        end = std::clock();
        time_combn[i] = static_cast<double>(end - begin);
      } else {
        time_combn[i] = std::numeric_limits<double>::infinity();
      }

      begin = std::clock();
      soln_glynn = permanent::glynn_square<T, I>(m, n, array);
      ensure(soln_glynn);
      end = std::clock();
      time_glynn[i] = static_cast<double>(end - begin);

      begin = std::clock();
      soln_ryser = permanent::ryser_square<T, I>(m, n, array);
      ensure(soln_ryser);
      end = std::clock();
      time_ryser[i] = static_cast<double>(end - begin);

      if (run_combinatorial) {
        if ((std::fabs(soln_combn - soln_glynn) /
                 std::min(std::fabs(soln_combn), std::fabs(soln_glynn)) >
             TOLERANCE) ||
            (std::fabs(soln_ryser - soln_combn) /
                 std::min(std::fabs(soln_ryser), std::fabs(soln_combn)) >
             TOLERANCE) ||
            (std::fabs(soln_glynn - soln_ryser) /
                 std::min(std::fabs(soln_glynn), std::fabs(soln_ryser)) >
             TOLERANCE)) {
          std::cerr << "Bad permanent values:"
                    << "\nCombn: " << soln_combn << "\nGlynn: " << soln_glynn
                    << "\nRyser: " << soln_ryser << std::endl;
          csv_file.close();
          return 3;
        }
      } else {
        // Only check Glynn vs Ryser when not running combinatorial
        if (std::fabs(soln_glynn - soln_ryser) /
                std::min(std::fabs(soln_glynn), std::fabs(soln_ryser)) >
            TOLERANCE) {
          std::cerr << "Bad permanent values:"
                    << "\nGlynn: " << soln_glynn
                    << "\nRyser: " << soln_ryser << std::endl;
          csv_file.close();
          return 3;
        }
      }
    }
  } else {
    for (size_t i = 0; i != NUM_TRIALS; ++i) {
      if (run_combinatorial) {
        begin = std::clock();
        soln_combn = permanent::combinatoric_rectangular<double>(m, n, array);
        ensure(soln_combn);
        end = std::clock();
        time_combn[i] = static_cast<double>(end - begin);
      } else {
        time_combn[i] = std::numeric_limits<double>::infinity();
      }

      begin = std::clock();
      soln_glynn = permanent::glynn_rectangular<double>(m, n, array);
      ensure(soln_glynn);
      end = std::clock();
      time_glynn[i] = static_cast<double>(end - begin);

      begin = std::clock();
      soln_ryser = permanent::ryser_rectangular<double>(m, n, array);
      ensure(soln_ryser);
      end = std::clock();
      time_ryser[i] = static_cast<double>(end - begin);

      if (run_combinatorial) {
        if ((std::fabs(soln_combn - soln_glynn) /
                 std::min(std::fabs(soln_combn), std::fabs(soln_glynn)) >
             TOLERANCE) ||
            (std::fabs(soln_ryser - soln_combn) /
                 std::min(std::fabs(soln_ryser), std::fabs(soln_combn)) >
             TOLERANCE) ||
            (std::fabs(soln_glynn - soln_ryser) /
                 std::min(std::fabs(soln_glynn), std::fabs(soln_ryser)) >
             TOLERANCE)) {
          std::cerr << std::scientific << "Bad permanent values:"
                    << "\nCombn: " << soln_combn << "\nGlynn: " << soln_glynn
                    << "\nRyser: " << soln_ryser << std::endl;
          csv_file.close();
          return 3;
        }
      } else {
        // Only check Glynn vs Ryser when not running combinatorial
        if (std::fabs(soln_glynn - soln_ryser) /
                std::min(std::fabs(soln_glynn), std::fabs(soln_ryser)) >
            TOLERANCE) {
          std::cerr << std::scientific << "Bad permanent values:"
                    << "\nGlynn: " << soln_glynn
                    << "\nRyser: " << soln_ryser << std::endl;
          csv_file.close();
          return 3;
        }
      }
    }
  }

  // Calculate the mean for the runtime of each algorithm
  double mean_combn = 0;
  double mean_glynn = 0;
  double mean_ryser = 0;
  for (size_t i = 0; i != NUM_TRIALS; ++i) {
    if (run_combinatorial) {
      mean_combn += time_combn[i];
    }
    mean_glynn += time_glynn[i];
    mean_ryser += time_ryser[i];
  }
  if (run_combinatorial) {
    mean_combn /= NUM_TRIALS;
  } else {
    mean_combn = std::numeric_limits<double>::infinity();
  }
  mean_glynn /= NUM_TRIALS;
  mean_ryser /= NUM_TRIALS;

  int fastest;
  if (mean_ryser <= mean_glynn) {
    fastest = 0;
  } else if (mean_combn <= mean_glynn) {
    fastest = 1;
  } else {
    fastest = 2;
  }

  // Write line
  std::cout << std::scientific << static_cast<double>(m) / n << ',' << std::setw(2) << n << ','
            << std::scientific << mean_combn << ',' << mean_glynn << ',' << std::scientific
            << mean_ryser << ',' << fastest << std::endl;

  csv_file << std::scientific << static_cast<double>(m) / n << ',' << std::setw(2) << n << ','
           << std::scientific << mean_combn << ',' << mean_glynn << ',' << std::scientific
           << mean_ryser << ',' << fastest << '\n';

  if (csv_file.fail()) {
    std::cerr << "Error writing to CSV file" << std::endl;
    csv_file.close();
    return 2;
  }
  return 0;
}

template <typename T, typename I = void>
int run_tuning(const char *filename)
{
  // Open CSV file
  std::ofstream csv_file = std::ofstream(filename);
  // csv_file.open(filename);
  if (csv_file.fail()) {
    std::cerr << "Cannot open CSV file `" << filename << '`' << std::endl;
    return 2;
  }

  // Print CSV headers
  csv_file << CSV_HEADER << '\n';
  if (csv_file.fail()) {
    std::cerr << "Error writing to CSV file `" << filename << '`' << std::endl;
    csv_file.close();
    return 2;
  }

  // Set format options for printing
  csv_file.precision(9);
  csv_file.fill(' ');

  // Random binary matrix for testing
  double array[MAX_ROWS * MAX_COLS];
  std::srand(static_cast<unsigned int>(0x23a23cf5033c3c81UL));
  for (size_t i = 0; i < MAX_ROWS * MAX_COLS; ++i) {
    array[i] = (std::rand() / static_cast<double>(RAND_MAX) - 0.5) * 2;
  }

  double time_combn[128];
  double time_glynn[128];
  double time_ryser[128];

  // Iterate over number of rows and number of columns
  for (size_t m = 2; m <= MAX_ROWS; ++m) {
    for (size_t n = m; n <= MAX_COLS; ++n) {
      int result =
          generate_tuning_data<T, I>(csv_file, m, n, array, time_combn, time_glynn, time_ryser);
      if (result != 0) {
        return result;
      }
    }
  }

  // Close CSV file
  csv_file.close();

  // Exit successfully
  return 0;
}

}  // namespace

}  // namespace permanent

int main(int argc, const char **argv)
{
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " CSV_FILE" << std::endl;
    return 1;
  }
  return permanent::run_tuning<double>(argv[1]);
}

#undef _PERMANENT_ALWAYS_INLINE
#undef _PERMANENT_DEFAULT_TUNING
