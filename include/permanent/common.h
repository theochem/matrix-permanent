/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_common_h_)
#define permanent_common_h_

#include <complex>
#include <cstdint>
#include <type_traits>

namespace permanent {

// Numerical types

using size_t = std::size_t;

using sgn_t = std::int_fast8_t;

// Static dispatcher for result type

template <typename T, typename = void>
struct _result_t;

template <typename T>
struct _result_t<T, std::enable_if_t<std::is_scalar_v<T>, void>>
{
  typedef double type;
};

template <typename T>
struct _result_t<std::complex<T>, std::enable_if_t<std::is_scalar_v<T>, void>>
{
  typedef std::complex<double> type;
};

template <typename T>
using result_t = _result_t<T>::type;

// Compute the parity of an integer

template <typename T>
inline auto parity(const T &x)
{
  return 1 - ((x & 1) << 1);
};

}  // namespace permanent

#endif  // permanent_common_h_
