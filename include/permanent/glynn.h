/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_glynn_h_)
#define permanent_glynn_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/tables.h>

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> glynn_square(const size_t m, const size_t, const T *ptr)
{
  using namespace complex_ops;

  sgn_t delta[64];        // delta array (all +1 to start)
  size_t perm[64 + 1];    // permutation array with {0, ..., m}
  result_t<T, I> vec[64]; // intermediate row sums
  result_t<T, I> out = 1; // result
  for (size_t i = 0; i != m; ++i) {
    delta[i] = 1;
  }
  for (size_t i = 0; i != m; ++i) {
    perm[i] = i;
  }

  // first permutation
  for (size_t j = 0; j != m; ++j) {
    result_t<T, I> sum = 0;
    for (size_t i = 0; i != m; ++i) {
      sum += ptr[i * m + j] * delta[i];
    }
    vec[j] = sum;
  }
  for (size_t j = 0; j != m; ++j) {
      out *= vec[j];
  }

  // further permutations
  size_t bound = m - 1;
  size_t pos = 0;
  sgn_t sign = 1;
  while (pos != bound) {
    // update sign and delta
    sign *= -1;
    delta[bound - pos] *= -1;
    // compute term
    for (size_t j = 0; j != m; ++j) {
      vec[j] += delta[bound - pos] * ptr[m * (bound - pos) + j] * 2;
    }
    result_t<T, I> prod = 1;
    for (size_t i = 0; i < m; ++i) {
        prod *= vec[i];
    }
    out += sign * prod;
    // advance permutation
    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
  }

  return out / (static_cast<size_t>(1) << bound);
}

template <typename T, typename I = void>
result_t<T, I> glynn_rectangular(const size_t m, const size_t n, const T *ptr)
{
  using namespace complex_ops;

  sgn_t delta[64];        // delta array (all +1 to start)
  size_t perm[64 + 1];    // permutation array with {0, ..., m}
  result_t<T, I> vec[64]; // intermediate row sums
  result_t<T, I> out = 1; // result
  for (size_t i = 0; i != n; ++i) {
    delta[i] = 1;
  }
  for (size_t i = 0; i != n; ++i) {
    perm[i] = i;
  }

  // first permutation
  for (size_t j = 0; j != n; ++j) {
    result_t<T, I> sum = 0;
    for (size_t i = 0; i != m; ++i) {
      sum += ptr[i * n + j] * delta[i];
    }
    for (size_t i = m; i != n; ++i) {
      sum += delta[i];
    }
    vec[j] = sum;
  }
  for (size_t j = 0; j != n; ++j) {
      out *= vec[j];
  }

  // further permutations
  size_t pos = 0;
  size_t bound = n - 1;
  sgn_t sign = 1;
  while (pos != bound) {
    // update sign and delta
    sign *= -1;
    delta[bound - pos] *= -1;
    // compute term
    if (bound - pos < m) {
      for (size_t j = 0; j != n; ++j) {
        vec[j] += delta[bound - pos] * ptr[n * (bound - pos) + j] * 2;
      }
    } else {
      for (size_t j = 0; j != n; ++j) {
        vec[j] += delta[bound - pos] * 2;
      }
    }
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != n; ++i) {
        prod *= vec[i];
    }
    out += sign * prod;

    // advance permutation
    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
  }

  return (out / (static_cast<size_t>(1) << bound)) / factorial<result_t<T, I>>(n - m);
}

template <typename T, typename I = void>
result_t<T, I> glynn(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? glynn_square<T, I>(m, n, ptr) : glynn_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_glynn_h_