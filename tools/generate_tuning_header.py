import sys
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.metrics import accuracy_score

CSV_FILE = sys.argv[1]
HEADER_FILE = sys.argv[2]

# Read output from tuning program
tuning = pd.read_csv(CSV_FILE, usecols=["M/N", "N", "Fastest"])

# Split data into two parts based on N=13
tuning_small = tuning[tuning["N"] <= 13].copy()
tuning_large = tuning[tuning["N"] > 13].copy()

# Update label columns for ML
tuning_small.rename(columns={"N": "x", "M/N": "y", "Fastest": "target"}, inplace=True)
tuning_large.rename(columns={"N": "x", "M/N": "y", "Fastest": "target"}, inplace=True)

# Find Combinatoric Square limit
matching_row = tuning["Fastest"] == 1
combn_limit = tuning.loc[matching_row, "N"].iloc[-1] if matching_row.any() else 0

# Update classes for first SVM (N <= 13, all three algorithms)
update_ryser = tuning_small["target"] == 0
tuning_small.loc[update_ryser, "target"] = -1
update_glynn = tuning_small["target"] == 2
tuning_small.loc[update_glynn, "target"] = -1

# Update classes for second SVM (N > 13, Glynn vs Ryser)
tuning_large["target"] = np.where(tuning_large["target"] == 0, -1, 1)

# Train first SVM (N <= 13)
features_small = tuning_small[["x", "y"]]
label_small = tuning_small["target"]
size_small = tuning_small.shape[0]
test_size_small = int(np.round(size_small * 0.1, 0))

x_train_small = features_small[:-test_size_small].values
y_train_small = label_small[:-test_size_small].values
x_test_small = features_small[-test_size_small:].values
y_test_small = label_small[-test_size_small:].values

linear_model_small = svm.SVC(kernel="linear", C=100.0)
linear_model_small.fit(x_train_small, y_train_small)

# Train second SVM (N > 13)
features_large = tuning_large[["x", "y"]]
label_large = tuning_large["target"]
size_large = tuning_large.shape[0]
test_size_large = int(np.round(size_large * 0.1, 0))

x_train_large = features_large[:-test_size_large].values
y_train_large = label_large[:-test_size_large].values
x_test_large = features_large[-test_size_large:].values
y_test_large = label_large[-test_size_large:].values

linear_model_large = svm.SVC(kernel="linear", C=100.0)
linear_model_large.fit(x_train_large, y_train_large)

# Get coefficients and bias for both models
coefficients_small = linear_model_small.coef_[0]
bias_small = linear_model_small.intercept_[0]
coefficients_large = linear_model_large.coef_[0]
bias_large = linear_model_large.intercept_[0]

# Write header file with all parameters
param_1 = coefficients_small[0]  # First hyperplane x coefficient
param_2 = coefficients_small[1]  # First hyperplane y coefficient
param_3 = bias_small            # First hyperplane bias
param_4 = combn_limit          # Combinatoric square limit
param_5 = coefficients_large[0] # Second hyperplane x coefficient
param_6 = coefficients_large[1] # Second hyperplane y coefficient
param_7 = bias_large           # Second hyperplane bias
param_8 = 13.0                 # Combinatorial limit

try:
    with open(HEADER_FILE, "w") as file_ptr:
        file_ptr.write(
            f"""/* Copyright 2024 QC-Devs (GPLv3) */

#if !defined(permanent_tuning_h_)
#define permanent_tuning_h_

namespace permanent {{

template<typename Type, typename IntType = void>
struct _tuning_params_t
{{
  static constexpr double PARAM_1 = {param_1:+16.9e};
  static constexpr double PARAM_2 = {param_2:+16.9e};
  static constexpr double PARAM_3 = {param_3:+16.9e};
  static constexpr double PARAM_4 = {param_4:+16.9e};
  static constexpr double PARAM_5 = {param_5:+16.9e};
  static constexpr double PARAM_6 = {param_6:+16.9e};
  static constexpr double PARAM_7 = {param_7:+16.9e};
  static constexpr double PARAM_8 = {param_8:+16.9e};
}};

template<typename Type, typename IntType = void>
static constexpr double PARAM_1 = _tuning_params_t<Type, IntType>::PARAM_1;

template<typename Type, typename IntType = void>
static constexpr double PARAM_2 = _tuning_params_t<Type, IntType>::PARAM_2;

template<typename Type, typename IntType = void>
static constexpr double PARAM_3 = _tuning_params_t<Type, IntType>::PARAM_3;

template<typename Type, typename IntType = void>
static constexpr double PARAM_4 = _tuning_params_t<Type, IntType>::PARAM_4;

template<typename Type, typename IntType = void>
static constexpr double PARAM_5 = _tuning_params_t<Type, IntType>::PARAM_5;

template<typename Type, typename IntType = void>
static constexpr double PARAM_6 = _tuning_params_t<Type, IntType>::PARAM_6;

template<typename Type, typename IntType = void>
static constexpr double PARAM_7 = _tuning_params_t<Type, IntType>::PARAM_7;

template<typename Type, typename IntType = void>
static constexpr double PARAM_8 = _tuning_params_t<Type, IntType>::PARAM_8;

}}  // namespace permanent

#endif  // permanent_tuning_h_
"""
        )

except IOError:
    print("Cannot open file!")
    exit(1)

except Exception as e:
    print("Error occurred:", e)
    exit(1)