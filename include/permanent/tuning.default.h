/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_tuning_h_)
#define permanent_tuning_h_

namespace permanent {

template <typename Type, typename IntType = void>
struct _tuning_params_t
{
  // First hyperplane parameters (N ≤ 13): -6.67(m/n) + (-0.35)(n) + 6.49 = 0
  static constexpr double PARAM_1 = -6.67000e+00;  // coefficient of m/n
  static constexpr double PARAM_2 = -3.50000e-01;  // coefficient of n
  static constexpr double PARAM_3 = +6.49000e+00;  // constant term
  static constexpr double PARAM_4 = +4.00000e+00;  // Combinatorial crossover limit

  // Second hyperplane parameters (N > 13): -10.76(m/n) + 0.09(n) + 1.96 = 0
  static constexpr double PARAM_5 = -1.07600e+01;  // coefficient of m/n
  static constexpr double PARAM_6 = +9.00000e-02;  // coefficient of n
  static constexpr double PARAM_7 = +1.96000e+00;  // constant term
  static constexpr double PARAM_8 = +1.30000e+01;  // Combinatorial limit (N=13)
};

template <typename Type, typename IntType = void>
static constexpr double PARAM_1 = _tuning_params_t<Type, IntType>::PARAM_1;

template <typename Type, typename IntType = void>
static constexpr double PARAM_2 = _tuning_params_t<Type, IntType>::PARAM_2;

template <typename Type, typename IntType = void>
static constexpr double PARAM_3 = _tuning_params_t<Type, IntType>::PARAM_3;

template <typename Type, typename IntType = void>
static constexpr double PARAM_4 = _tuning_params_t<Type, IntType>::PARAM_4;

template <typename Type, typename IntType = void>
static constexpr double PARAM_5 = _tuning_params_t<Type, IntType>::PARAM_5;

template <typename Type, typename IntType = void>
static constexpr double PARAM_6 = _tuning_params_t<Type, IntType>::PARAM_6;

template <typename Type, typename IntType = void>
static constexpr double PARAM_7 = _tuning_params_t<Type, IntType>::PARAM_7;

template <typename Type, typename IntType = void>
static constexpr double PARAM_8 = _tuning_params_t<Type, IntType>::PARAM_8;

}  // namespace permanent

#endif  // permanent_tuning_h_