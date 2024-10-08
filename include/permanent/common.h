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

template <typename Type, typename IntType = void, typename = void>
struct _result_t;

template <typename Type>
struct _result_t<Type, void, std::enable_if_t<std::is_scalar_v<Type>, void>>
{
  typedef double type;
};

template <typename Type>
struct _result_t<std::complex<Type>, void, std::enable_if_t<std::is_scalar_v<Type>, void>>
{
  typedef std::complex<double> type;
};

template <typename Type, typename IntType>
struct _result_t<Type, IntType,
                 std::enable_if_t<std::is_integral_v<Type> && std::is_integral_v<IntType>, void>>
{
  typedef IntType type;
};

template <typename Type, typename IntType>
struct _result_t<std::complex<Type>, std::complex<IntType>,
                 std::enable_if_t<std::is_integral_v<Type> && std::is_integral_v<IntType>, void>>
{
  typedef std::complex<IntType> type;
};

template <typename Type, typename IntType = void>
using result_t = _result_t<Type, IntType>::type;

}  // namespace permanent

#endif  // permanent_common_h_
