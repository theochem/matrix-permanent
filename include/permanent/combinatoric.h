/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_combinatoric_h_)
#define permanent_combinatoric_h_

#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/kperm-gray.h>
#include <permanent/perm-mv0.h>

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> combinatoric_square(const size_t m, const size_t n, const T *ptr)
{
  (void)n;

  perm_mv0 permutations(m);

  result_t<T, I> out = 0;
  const size_t *perm = permutations.data();
  do {
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != m; ++i) {
      prod *= ptr[i * m + perm[i]];
    }
    out += prod;
  } while (permutations.next());

  return out;
}

template <typename T, typename I = void>
result_t<T, I> combinatoric_rectangular(const size_t m, const size_t n, const T *ptr)
{
  kperm_gray permutations(n);
  permutations.first(m);

  result_t<T, I> out = 0;
  const size_t *perm = permutations.data();
  do {
    result_t<T, I> prod = 1;
    for (size_t i = 0; i != m; ++i) {
      prod *= ptr[i * n + perm[i]];
    }
    out += prod;
  } while (permutations.next());

  return out;
}

template <typename T, typename I = void>
result_t<T, I> combinatoric(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? combinatoric_square<T, I>(m, n, ptr)
                  : combinatoric_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_combinatoric_h_
