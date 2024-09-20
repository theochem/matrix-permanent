/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(PERMANENT_COMMON_H_)
#define PERMANENT_COMMON_H_

#include <cstdint>

#include <complex>
#include <type_traits>

namespace permanent {

// Tuning file

#if defined(WITH_TUNING_FILE)
// Include generated tuning parameters
#include "tuning.h"
#else
// Set default tuning parameters
constexpr double PARAM_1 = -0.572098;
constexpr double PARAM_2 = -22.014212;
constexpr double PARAM_3 = 15.29794;
constexpr double PARAM_4 = 3.0;
#endif

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
inline auto parity(const T &x) {
    return 1 - ((x & 1) << 1);
};

}  // namespace permanent

#endif  // PERMANENT_COMMON_H_
