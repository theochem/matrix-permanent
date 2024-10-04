/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(PERMANENT_COMPLEX_H_)
#define PERMANENT_COMPLEX_H_

#include <complex>

namespace permanent {

// Allow type promotion to complex

template <typename T>
struct _identity_t
{
  typedef T type;
};

template <typename T>
using identity_t = _identity_t<T>::type;

}  // namespace permanent

#define DEF_COMPLEX_OPS(OP)                                                                \
                                                                                           \
  template <typename T>                                                                    \
  std::complex<T> operator OP(const std::complex<T> &x, const permanent::identity_t<T> &y) \
  {                                                                                        \
    return x OP y;                                                                         \
  }                                                                                        \
                                                                                           \
  template <typename T>                                                                    \
  std::complex<T> operator OP(const permanent::identity_t<T> &x, const std::complex<T> &y) \
  {                                                                                        \
    return x OP y;                                                                         \
  }

DEF_COMPLEX_OPS(+)
DEF_COMPLEX_OPS(-)
DEF_COMPLEX_OPS(*)
DEF_COMPLEX_OPS(/)

#undef DEF_COMPLEX_OPS

#endif  // PERMANENT_COMPLEX_H_
