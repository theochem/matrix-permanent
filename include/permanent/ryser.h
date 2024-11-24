/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_ryser_h_)
#define permanent_ryser_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/tables.h>

#include <bit>
#include <vector>

namespace {

inline void all_combinations(int n, int k, std::vector<std::vector<int>> &k_perms)
{
  k_perms.clear();  // Clear k_perms to make sure it starts empty
  k_perms.reserve(1 << n);

  for (int i = 0; i < (1 << n); i++) {
    unsigned int cur = i ^ (i >> 1);
    if (std::popcount(cur) == k) {
      std::vector<int> vec;  // Initialize vector for new combination
      vec.reserve(n);
      for (int j = 0; j < n; j++) {
        if (cur & (1 << j)) {
          vec.push_back(j);  // Add element to vector
        }
      }
      k_perms.push_back(vec);  // Store the combination
    }
  }
}

template <typename T>
auto parity(const T &x)
{
  return 1 - ((x & 1) << 1);
};

template <typename T>
void swap2(T *perm, T i, T j)
{
  const T temp = perm[i];
  perm[i] = perm[j];
  perm[j] = temp;
}

template <typename T>
void init_perm(const T N, T *const the_fact_set, T *const the_perm_set, T *const the_inv_set)
{
  for (T i = 0; i < N; i++) {
    the_fact_set[i + 1] = 0;
    the_perm_set[i] = i;
    the_inv_set[i] = i;
  }
}

template <typename T>
bool gen_next_perm(T *const falling_fact, T *const perm, T *const inv_perm, const T cols,
                   const T u_)
{
  /* Use the current permutation to check for next permutation
   * lexicographically, if it exists update the curr_array by swapping the
   * leftmost changed element (at position i < k). Replace the elements up to
   * position k by the smallest element lying to the right of i. Do not put
   * remaining k, .., n - 1 positions in ascending order to improve efficiency
   * has time complexity O(n) in the worst case. */
  T i = u_ - 1;
  T m1 = cols - i - 1;

  /* Begin update of falling_fact - recall falling_fact[0] = 0, so increment
   * index accordingly. If i becomes negative during the check, you are
   * pointing to the sentinel value, so you are done generating
   * permutations. */
  while (falling_fact[i + 1] == m1) {
    falling_fact[i + 1] = 0;
    ++m1;
    --i;
  }
  if (i == 0UL - 1) {
    return false;
  }
  ++falling_fact[i + 1];

  /* Find smallest element perm[i] < perm[j] that lies to the right of
   * pos i, and then update the state of the permutation using its inverse
   * to generate next. */
  T z = perm[i];
  do {
    ++z;
  } while (inv_perm[z] <= i);
  const T j = inv_perm[z];
  swap2(perm, i, j);
  swap2(inv_perm, perm[i], perm[j]);
  ++i;
  T y = 0;

  /* Find the smallest elements to the right of position i. */
  while (i < u_) {
    while (inv_perm[y] < i) {
      ++y;
    }
    const T k = inv_perm[y];
    swap2(perm, i, k);
    swap2(inv_perm, perm[i], perm[k]);
    ++i;
  }
  return true;
}

}  // namespace

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> ryser_square(const size_t m, const size_t, const T *ptr)
{
  T out = 0;
  size_t c = 1UL << m;

  for (size_t k = 0; k < c; ++k) {
    T rowsumprod = 1.0;
    for (size_t i = 0; i < m; ++i) {
      T rowsum = 0.0;
      for (size_t j = 0; j < m; ++j) {
        if (k & (1UL << j)) {
          rowsum += ptr[m * i + j];
        }
      }
      rowsumprod *= rowsum;
    }
    out += rowsumprod * (1 - ((std::popcount(k) & 1) << 1));
  }

  return static_cast<result_t<T, I>>(out * ((m % 2 == 1) ? -1 : 1));
}

template <typename T, typename I = void>
result_t<T, I> ryser_rectangular(const size_t m, const size_t n, const T *ptr)
{
  int sign = 1;
  result_t<T, I> out = 0.0;
  std::vector<std::vector<int>> k_perms;

  /* Iterate over subsets from size 0 to size m */

  for (size_t k = 0; k < m; ++k) {
    all_combinations(n, m - k, k_perms);
    size_t bin = binomial((n - m + k), k);
    result_t<T, I> permsum = 0.0;

    /* Compute permanents of each submatrix */
    result_t<T, I> vec[64];
    for (const auto &combination : k_perms) {
      result_t<T, I> colprod;
      for (size_t i = 0; i < m; ++i) {
        result_t<T, I> matsum = 0.0;
        for (size_t j = 0; j < (m - k); ++j) matsum += ptr[i * n + combination[j]];
        vec[i] = matsum;
      }

      colprod = 1.0;
      for (size_t i = 0; i < m; ++i) colprod *= vec[i];
      permsum += colprod * sign * bin;
    }

    /* Add term to result */

    out += permsum;

    sign *= -1;
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
