/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_opt_h_)
#define permanent_opt_h_

#include <permanent/combinatoric.h>
#include <permanent/common.h>
#include <permanent/complex.h>
#include <permanent/glynn.h>
#include <permanent/ryser.h>

#if defined(_PERMANENT_DEFAULT_TUNING)
#include <permanent/tuning.default.h>
#else
#include <permanent/tuning.h>
#endif

namespace permanent {

template <typename T, typename I = void>
result_t<T, I> opt_square(const size_t m, const size_t n, const T *ptr)
{
  if (n <= PARAM_4<T, I>) {
    return ryser_square<T, I>(m, n, ptr);
  } else if (n * PARAM_1<T, I> + PARAM_2<T, I> + PARAM_3<T, I> > 0) {
    return combinatoric_square<T, I>(m, n, ptr);
  } else {
    return glynn_square<T, I>(m, n, ptr);
  }
}

template <typename T, typename I = void>
result_t<T, I> opt_rectangular(const size_t m, const size_t n, const T *ptr)
{
  if (n * PARAM_1<T, I> + PARAM_2<T, I> * m / n + PARAM_3<T, I> > 0) {
    return combinatoric_rectangular<T, I>(m, n, ptr);
  } else {
    return glynn_rectangular<T, I>(m, n, ptr);
  }
}

template <typename T, typename I = void>
result_t<T, I> opt(const size_t m, const size_t n, const T *ptr)
{
  return (m == n) ? opt_square<T, I>(m, n, ptr) : opt_rectangular<T, I>(m, n, ptr);
}

}  // namespace permanent

#endif  // permanent_opt_h_
