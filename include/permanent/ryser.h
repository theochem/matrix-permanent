/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_ryser_h_)
#define permanent_ryser_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/tables.h>

#include <bit>

namespace {

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
result_t<T, I> ryser_square(const size_t m, const size_t n, const T *ptr)
{
  (void)n;

  size_t i, j;
  size_t k;
  result_t<T, I> rowsum;
  result_t<T, I> out = 0;

  /* Iterate over c = pow(2, m) submatrices (equal to (1 << m)) submatrices.
   */

  size_t c = 1UL << m;

  /* Loop over columns of submatrix; compute product of row sums. */

  for (k = 0; k < c; ++k) {
    result_t<T, I> rowsumprod = 1.0;
    for (i = 0; i < m; ++i) {
      /* Loop over rows of submatrix; compute row sum. */

      rowsum = 0.0;
      for (j = 0; j < m; ++j) {
        /* Add element to row sum if the row index is in the
         * characteristic * vector of the submatrix, which is the binary
         * vector given by k.  */

        if (k & (1UL << j)) {
          rowsum += ptr[m * i + j];
        }
      }

      /* Update product of row sums. */

      rowsumprod *= rowsum;
    }

    /* Add term multiplied by the parity of the characteristic vector. */

    out += rowsumprod * parity(std::popcount(k));
  }

  /* Return answer with the correct sign (times -1 for odd m). */

  return out * parity(m);
}

template <typename T, typename I = void>
result_t<T, I> ryser_rectangular(const size_t m, const size_t n, const T *ptr)
{
  size_t falling_fact[128];
  size_t perm[128];
  size_t inv_perm[128];
  falling_fact[0] = 0;

  size_t i, j, k;

  int value_sign = 1;

  result_t<T, I> sum_of_matrix_vals;
  result_t<T, I> sum_over_k_vals = 0.0;
  result_t<T, I> vec[64];

  /* Dealing with a rectangle. Can't use bit hacking trick here. */
  for (k = 0; k < m; ++k) {
    /* Store the binomial coefficient for this k value bin_c. */
    size_t counter = 0;  // Count how many permutations you have generated
    size_t bin_c = binomial((n - m + k), k);
    result_t<T, I> prod_of_cols = 1;
    result_t<T, I> result = 0;

    /* (Re)initialize the set to permute for this k value. */
    init_perm(n, falling_fact, perm, inv_perm);

    /* sort up to position u + 1 where u = min(m - k, n - 1). */
    size_t sort_up_to = n - 1;

    if ((m - k) < sort_up_to) {
      sort_up_to = (m - k);
    }

    for (i = 0; i < m; ++i) {
      sum_of_matrix_vals = 0.0;
      for (j = 0; j < sort_up_to; ++j) {
        sum_of_matrix_vals += ptr[i * n + perm[j]];
      }
      vec[i] = sum_of_matrix_vals;
    }
    for (i = 0; i < m; ++i) {
      prod_of_cols *= vec[i];
    }

    result += value_sign * bin_c * prod_of_cols;
    counter += 1;

    /* Iterate over second to last permutations of the set. */
    while (gen_next_perm(falling_fact, perm, inv_perm, n, sort_up_to)) {
      for (i = 0; i < m; ++i) {
        sum_of_matrix_vals = 0.0;
        for (j = 0; j < m - k; ++j) {
          sum_of_matrix_vals += ptr[i * n + perm[j]];
        }
        vec[i] = sum_of_matrix_vals;
      }
      prod_of_cols = 1.0;
      for (i = 0; i < m; ++i) {
        prod_of_cols *= vec[i];
      }

      result += value_sign * bin_c * prod_of_cols;
      counter += 1;
    }
    sum_over_k_vals += (result / counter) * binomial(n, (m - k));
    value_sign *= -1;
  }
  return sum_over_k_vals;
}

template <typename T, typename I = void>
result_t<T, I> ryser(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? ryser_square<T, I>(m, n, ptr) : ryser_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_ryser_h_
