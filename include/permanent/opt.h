/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_opt_h_)
#define permanent_opt_h_

#include <permanent/common.h>
#include <permanent/complex.h>

#include <permanent/combinatoric.h>
#include <permanent/glynn.h>
#include <permanent/ryser.h>

#if defined(_PERMANENT_DEFAULT_TUNING)
#include <permanent/tuning.default.h>
#else
#include <permanent/tuning.h>
#endif

namespace permanent {

// Polynomial feature computation helper
// Generates features in sklearn PolynomialFeatures order:
// For degree d: N^d, N^(d-1)*(M/N), N^(d-2)*(M/N)^2, ..., (M/N)^d
template <typename Type>
inline void compute_poly_features(const size_t n, const double ratio, double* features)
{
  int idx = 0;

  // Generate all combinations: N^i * (M/N)^j where i+j <= POLY_DEGREE
  // sklearn uses lexicographic order: higher powers of first feature come first
  for (int total_degree = 0; total_degree <= POLY_DEGREE; ++total_degree) {
    for (int n_power = total_degree; n_power >= 0; --n_power) {
      int mn_power = total_degree - n_power;

      double term = 1.0;

      // Compute N^n_power
      for (int i = 0; i < n_power; ++i) {
        term *= static_cast<double>(n);
      }

      // Compute (M/N)^mn_power
      for (int i = 0; i < mn_power; ++i) {
        term *= ratio;
      }

      features[idx++] = term;
    }
  }
}

// Predict optimal algorithm using polynomial logistic regression
// Returns: 0 = Ryser, 1 = Glynn
template <typename Type>
inline int predict_algorithm(const size_t m, const size_t n)
{
  const double ratio = static_cast<double>(m) / static_cast<double>(n);

  // 1. Compute polynomial features
  double features[N_FEATURES];
  compute_poly_features<Type>(n, ratio, features);

  // 2. Scale features
  double scaled_features[N_FEATURES];
  for (int i = 0; i < N_FEATURES; ++i) {
    scaled_features[i] = (features[i] - SCALER_MEAN[i]) / SCALER_SCALE[i];
  }

  // 3. Compute decision values for each class
  double decision_values[N_CLASSES];
  for (int c = 0; c < N_CLASSES; ++c) {
    decision_values[c] = INTERCEPTS[c];
    for (int f = 0; f < N_FEATURES; ++f) {
      decision_values[c] += COEFFICIENTS[c * N_FEATURES + f] * scaled_features[f];
    }
  }

  // 4. Return class with maximum decision value
  int max_idx = 0;
  double max_val = decision_values[0];
  for (int i = 1; i < N_CLASSES; ++i) {
    if (decision_values[i] > max_val) {
      max_val = decision_values[i];
      max_idx = i;
    }
  }

  return max_idx;
}

template <typename T, typename I = void>
result_t<T, I> opt_square(const size_t m, const size_t n, const T *ptr)
{
  const int algo = predict_algorithm<T>(m, n);

  switch (algo) {
    case 0:  // Ryser
      return ryser_square<T, I>(m, n, ptr);
    case 1:  // Glynn
      return glynn_square<T, I>(m, n, ptr);
    default:
      // Fallback to Glynn for safety
      return glynn_square<T, I>(m, n, ptr);
  }
}

template <typename T, typename I = void>
result_t<T, I> opt_rectangular(const size_t m, const size_t n, const T *ptr)
{
  const int algo = predict_algorithm<T>(m, n);

  switch (algo) {
    case 0:  // Ryser
      return ryser_rectangular<T, I>(m, n, ptr);
    case 1:  // Glynn
      return glynn_rectangular<T, I>(m, n, ptr);
    default:
      // Fallback to Glynn for safety
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
