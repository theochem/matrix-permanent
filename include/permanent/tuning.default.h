/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_tuning_h_)
#define permanent_tuning_h_

namespace permanent {

template <typename Type, typename IntType = void>
struct _tuning_params_t
{
  static constexpr double PARAM_1 = -5.72098e-01;
  static constexpr double PARAM_2 = -2.20142e+01;
  static constexpr double PARAM_3 = +1.52979e+01;
  static constexpr double PARAM_4 = +3.00000e+00;
};

template <typename Type, typename IntType = void>
static constexpr double PARAM_1 = _tuning_params_t<Type, IntType>::PARAM_1;

template <typename Type, typename IntType = void>
static constexpr double PARAM_2 = _tuning_params_t<Type, IntType>::PARAM_2;

template <typename Type, typename IntType = void>
static constexpr double PARAM_3 = _tuning_params_t<Type, IntType>::PARAM_3;

template <typename Type, typename IntType = void>
static constexpr double PARAM_4 = _tuning_params_t<Type, IntType>::PARAM_4;

}  // namespace permanent

#endif  // permanent_tuning_h_
