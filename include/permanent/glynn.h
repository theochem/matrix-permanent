/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_glynn_h_)
#define permanent_glynn_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/tables.h>

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> glynn_square(const size_t m, const size_t n, const T *ptr)
{
  (void)n;

  /* Fill delta array (all +1 to start), and permutation array with [0...m]. */

  sgn_t delta[64];
  size_t perm[64 + 1];
  for (size_t i = 0; i < m; ++i) {
    perm[i] = i;
    delta[i] = 1;
  }

  /* Handle first permutation. */

  // result_t<T, I> vec[64];
  result_t<T, I> out = 1;
  for (size_t j = 0; j < m; ++j) {
    result_t<T, I> sum = 0;
    for (size_t i = 0; i < m; ++i) {
      sum += ptr[i * m + j] * delta[i];
    }
    out *= sum;
    // vec[j] = sum;
  }
  // for (size_t j = 0; j < m; ++j) {
  //     out *= vec[j];
  // }

  /* Iterate from the second to the final permutation. */

  std::size_t bound = m - 1;
  std::size_t pos = 0;
  sgn_t sign = 1;
  while (pos != bound) {
    /* Update sign and delta. */

    sign *= -1;
    delta[bound - pos] *= -1;

    /* Compute term. */

    result_t<T, I> prod = 1;
    for (size_t j = 0; j < m; ++j) {
      result_t<T, I> sum = 0;
      for (size_t i = 0; i < m; ++i) {
        sum += ptr[i * m + j] * delta[i];
      }
      // vec[j] = sum;
      prod *= sum;
    }

    /* Multiply by the product of the vectors in delta. */

    // for (size_t i = 0; i < m; ++i) {
    //     prod *= vec[i];
    // }
    out += sign * prod;

    /* Go to next permutation. */

    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
  }

  /* Divide by external factor and return permanent. */

  return out / (1UL << bound);
}

template <typename T, typename I = void>
result_t<T, I> glynn_rectangular(const size_t m, const size_t n, const T *ptr)
{
  /* Fill delta array (all +1 to start), and permutation array with [0...n].
   */

  sgn_t delta[64];
  size_t perm[64 + 1];
  for (size_t i = 0; i < n; ++i) {
    delta[i] = 1;
    perm[i] = i;
  }

  /* Handle first permutation. */

  result_t<T, I> out = 1;
  for (size_t j = 0; j < n; ++j) {
    result_t<T, I> sum = 0;
    for (size_t i = 0; i < m; ++i) {
      sum += ptr[i * n + j] * delta[i];
    }
    for (size_t k = m; k < n; ++k) {
      sum += delta[k];
    }
    // vec[j] = sum;
    out *= sum;
  }
  // for (i = 0; i < n; ++i) {
  //     out *= vec[i];
  // }

  /* Iterate from the second to the final permutation. */

  size_t pos = 0;
  size_t bound = n - 1;
  sgn_t sign = 1;
  while (pos != bound) {
    /* Update sign and delta. */

    sign *= -1;
    delta[bound - pos] *= -1;

    /* Compute term. */

    result_t<T, I> prod = 1;
    for (size_t j = 0; j < n; ++j) {
      result_t<T, I> sum = 0;
      for (size_t i = 0; i < m; ++i) {
        sum += ptr[i * n + j] * delta[i];
      }

      for (size_t k = m; k < n; ++k) {
        sum += delta[k];
      }
      prod *= sum;
    }

    /* Multiply by the product of the vectors in delta. */

    // T prod = 1.0;
    // for (i = 0; i < n; ++i) {
    //     prod *= vec[i];
    // }
    out += sign * prod;

    /* Go to next permutation. */

    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
  }

  /* Divide by external factor and return permanent. */

  return (out / (1UL << bound)) / factorial<result_t<T, I>>(n - m);
}

template <typename T, typename I = void>
result_t<T, I> glynn_square1(const size_t m, const size_t n, const T *ptr)
{
  (void)n;

  // Fill delta array ([+1...+1]) and permutation array ([0...m])
  size_t perm[64 + 1];
  sgn_t delta[64];
  for (size_t i = 0; i != m; ++i) {
    perm[i] = i;
    delta[i] = 1;
  }

  // Compute first permutation
  result_t<T, I> vec[64];
  result_t<T, I> out = 1;
  for (size_t j = 0; j != m; ++j) {
    // Compute inner terms
    result_t<T, I> sum = 0;
    for (size_t i = 0; i != m; ++i) {
      sum += ptr[i * m + j];
    }
    vec[j] = sum;

    // Compute first outer term
    out *= sum;
  }

  // Iterate from the second to the final permutation
  size_t bound = m - 1;
  size_t pos = 0;
  sgn_t sign = -1;
  while (pos != bound) {
    // Update delta
    size_t idx = bound - pos;
    delta[idx] *= -1;

    // Update inner terms
    for (size_t i = 0; i != m; ++i) {
      vec[i] += ptr[idx * m + i] * delta[idx] * 2;
    }

    // Add product of inner terms to outer term
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != m; ++i) {
      prod *= vec[i];
    }
    out += sign * prod;

    // Go to next permutation
    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
    sign *= -1;
  }

  // Divide outer term by external factor and return permanent
  return out / (1UL << bound);
}

template <typename T, typename I = void>
result_t<T, I> glynn_rectangular1(const size_t m, const size_t n, const T *ptr)
{
  // Fill delta array ([+1...+1]) and permutation array ([0...m])
  size_t perm[64 + 1];
  sgn_t delta[64];
  for (size_t i = 0; i != n; ++i) {
    perm[i] = i;
    delta[i] = 1;
  }

  // Compute first permutation
  result_t<T, I> vec[64];
  result_t<T, I> out = n * (n - m);
  for (size_t j = 0; j != m; ++j) {
    // Compute inner terms
    result_t<T, I> sum = 0;
    for (size_t i = 0; i != n; ++i) {
      sum += ptr[j * n + i];
    }
    vec[j] = sum;

    // Compute first outer term
    out *= sum;
  }
  for (size_t j = m; j != n; ++j) {
    vec[j] = n;
  }

  // Iterate from the second to the final permutation
  size_t bound = n - 1;
  size_t pos = 0;
  sgn_t sign = -1;
  while (pos != bound) {
    // Update delta
    size_t idx = bound - pos;
    delta[idx] *= -1;

    // Update inner terms
    for (size_t i = 0; i != m; ++i) {
      vec[i] += ptr[i * n + idx] * delta[idx] * 2;
    }
    for (size_t i = m; i != n; ++i) {
      vec[i] += delta[idx] * 2;
    }

    // Add product of inner terms to outer term
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != n; ++i) {
      prod *= vec[i];
    }
    out += sign * prod;

    // Go to next permutation
    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
    sign *= -1;
  }

  // Divide outer term by external factor and return permanent
  return (out / (1UL << bound)) / factorial<result_t<T, I>>(n - m);
}

template <typename T, typename I = void>
result_t<T, I> glynn_rectangular2(const size_t m, const size_t n, const T *ptr)
{
  // Fill delta array ([+1...+1]) and permutation array ([0...m])
  size_t perm[64 + 1];
  sgn_t delta[64];
  for (size_t i = 0; i != n; ++i) {
    perm[i] = i;
    delta[i] = 1;
  }

  // Compute first permutation
  result_t<T, I> vec[64];
  result_t<T, I> out = 1;
  for (size_t j = 0; j != n; ++j) {
    // Compute inner terms
    result_t<T, I> sum = n - m;
    for (size_t i = 0; i != m; ++i) {
      sum += ptr[i * n + j];
    }
    vec[j] = sum;

    // Compute first outer term
    out *= sum;
  }

  // Iterate from the second to the final permutation
  size_t bound = n - 1;
  size_t pos = 0;
  sgn_t sign = -1;
  while (pos != bound) {
    // Update delta
    size_t idx = bound - pos;
    delta[idx] *= -1;

    // Update inner terms
    if (idx < m) {
      for (size_t i = 0; i != n; ++i) {
        vec[i] += ptr[idx * n + i] * delta[idx] * 2;
      }
    } else {
      for (size_t i = 0; i != n; ++i) {
        vec[i] += delta[idx] * 2;
      }
    }

    // Add product of inner terms to outer term
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != n; ++i) {
      prod *= vec[i];
    }
    out += sign * prod;

    // Go to next permutation
    perm[0] = 0;
    perm[pos] = perm[pos + 1];
    ++pos;
    perm[pos] = pos;
    pos = perm[0];
    sign *= -1;
  }

  // Divide outer term by external factor and return permanent
  return (out / (1UL << bound)) / factorial<result_t<T, I>>(n - m);
}

template <typename T, typename I = void>
result_t<T, I> glynn(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? glynn_square<T, I>(m, n, ptr) : glynn_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_glynn_h_
