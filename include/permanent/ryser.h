/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_ryser_h_)
#define permanent_ryser_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/tables.h>

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> ryser_square(const size_t m, const size_t, const T *ptr)
{
  using namespace complex_ops;

  size_t subset[65] = {};          // k-subset with {1,..,n}
  result_t<T, I> rowsums[64] = {}; // intermediate row sums
  result_t<T, I> prods[64] = {};   // intermediate products of row sums

  size_t k = 0;
  do {
    result_t<T, I> prod = 1;
    // update subset
    if (k & 1) { // odd
      if (subset[k] - 1 == subset[k - 1]) { // remove subset[k] - 1 from position k - 1
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] -= ptr[m * row + subset[k] - 1 - 1]);
        }
        subset[k - 1] = subset[k];
        --k;
      } else { // insert subset[k] - 1 as second last element
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] += ptr[m * row + subset[k] - 1 - 1]);
        }
        subset[k + 1] = subset[k];
        --subset[k];
        ++k;
      }
    } else { // even
      if (subset[k] == m) { // remove m from end
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] -= ptr[m * row + m - 1]);
        }
        --k;
      } else { // append m
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] += ptr[m * row + m - 1]);
        }
        ++k;
        subset[k] = m;
      }
    }
    prods[k - 1] += prod;
  } while (k);

  result_t<T, I> out = 0;
  for (size_t i = 0; i != m; ++i) {
    out += prods[i] * (1 - (((ptrdiff_t)(i + 1) & 1) << 1));
  }

  return out * (1 - (((ptrdiff_t)m & 1) << 1));//(m & 1 ? -1 : 1);
}

template <typename T, typename I = void>
result_t<T, I> ryser_rectangular(const size_t m, const size_t n, const T *ptr)
{
  using namespace complex_ops;

  size_t subset[65] = { 0, 1 }; // k-subset with {1,...,n}
  result_t<T, I> rowsums[64];           // intermediate row sums
  result_t<T, I> prods[64] = { 1 };     // intermediate products of row sums

  for (size_t row = 0; row != m; ++row) {
    prods[0] *= (rowsums[row] = ptr[n * row]);
  }

  for (size_t k = 1; subset[1] != n;) {
    result_t<T, I> prod = 1;
    if (k & 1) { // odd
      if (subset[k] == n) { // remove n from end
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] -= ptr[n * row + n - 1]);
        }
        --k;
      } else {
        if (k < m) { // append n
          for (size_t row = 0; row != m; ++row) {
            prod *= (rowsums[row] += ptr[n * row + n - 1]);
          }
          ++k;
          subset[k] = n;
        } else { // increment subset[j] -> subset[j] + 1
          for (size_t row = 0; row != m; ++row) {
            prod *= (rowsums[row] += ptr[n * row + subset[k] + 1 - 1] - ptr[n * row + subset[k] - 1]);
          }
          ++subset[k];
        }
      }
    } else { // even
      if (subset[k - 1] == subset[k] - 1) { // remove subset[j] - 1 from position j - 1
        for (size_t row = 0; row != m; ++row) {
          prod *= (rowsums[row] -= ptr[n * row + subset[k] - 1 - 1]);
        }
        subset[k - 1] = subset[k];
        --k;
      } else {
        --subset[k];
        if (k < m) { // insert subset[j] - 1 as second last element
          for (size_t row = 0; row != m; ++row) {
            prod *= (rowsums[row] += ptr[n * row + subset[k] - 1 - 0]);
          }
          subset[k + 1] = subset[k] + 1;
          ++k;
        } else { // decrement subset[j] -> sunset[j] - 1
          for (size_t row = 0; row != m; ++row) {
            prod *= (rowsums[row] += ptr[n * row + subset[k] - 1 - 0] - ptr[n * row + subset[k] - 0]);
          }
        }
      }
    }
    prods[k - 1] += prod;
  }

  result_t<T, I> out = 0;
  for (size_t i = 0; i != m; ++i) {
    out += (prods[i] * binomial(n - i - 1, m - i - 1)) * (1 - (((ptrdiff_t)(m - i - 1) & 1) << 1));
  }

  return out;
}

template <typename T, typename I = void>
result_t<T, I> ryser(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? ryser_square<T, I>(m, n, ptr) : ryser_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_ryser_h_