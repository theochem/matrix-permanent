import numpy as np

import pandas as pd

from sklearn import svm
from sklearn.metrics import accuracy_score


CSV_FILE = "src/tuning.csv"

HEADER_FILE = "src/tuning.h"


# Read output from tuning program
tuning = pd.read_csv(CSV_FILE, usecols=["M/N", "N", "Fastest"])

# Update label columns for ML
tuning.rename(columns={"N": "x", "M/N": "y", "Fastest": "target"}, inplace=True)

# Find Ryser limit and update to dual class for SVM
# Locate rows where column value matches specified value
matching_row = tuning.loc[tuning["target"] == 0, "x"]

if not matching_row.empty:
    # Pull out specific column value from the last matching row
    ryser_limit = matching_row.iloc[-1].copy()
else:
    ryser_limit = 0

update_target = tuning["target"] == 0
tuning.loc[update_target, "target"] = -1

# Update classes to -1/1, Combn = 1
update_glynn = tuning["target"] == 2
tuning.loc[update_glynn, "target"] = -1

features = tuning[["x", "y"]]
label = tuning["target"]
value_counts = tuning["target"].value_counts()

# Create train/test split
size = tuning.shape[0]
test_size = int(np.round(size * 0.1, 0))

# Split dataset into training and testing sets
x_train = features[:-test_size].values
y_train = label[:-test_size].values
x_test = features[-test_size:].values
y_test = label[-test_size:].values

### Train a linear model - add a hard margin
linear_model = svm.SVC(kernel="linear", C=100.0)
linear_model.fit(x_train, y_train)

# Create grid to evaluate model
xx = np.linspace(-1, max(features["x"]) + 1, len(x_train))
yy = np.linspace(0, max(features["y"]) + 1, len(y_train))
YY, XX = np.meshgrid(yy, xx)
xy = np.vstack([XX.ravel(), YY.ravel()]).T
train_size = len(features[:-test_size]["x"])

# Get the separating hyperplane
Z = linear_model.decision_function(xy).reshape(XX.shape)

# Check the accuracy of the model
predictions_linear = linear_model.predict(x_test)
accuracy_linear = accuracy_score(y_test, predictions_linear)

# Get the coefficients and bias
coefficients = linear_model.coef_[0]
bias = linear_model.intercept_[0]

# Write a header file with constants defined as macros
param_1 = coefficients[0]
param_2 = coefficients[1]
param_3 = bias
param_4 = ryser_limit

try:
    with open(HEADER_FILE, "w") as file_ptr:
        file_ptr.write(
            f"""#ifndef TUNING_H_
#define TUNING_H_


 constexpr double PARAM_1 = {param_1:.9f};
 constexpr double PARAM_2 = {param_2:.9f};
 constexpr double PARAM_3 = {param_3:.9f};
 constexpr double PARAM_4 = {param_4:.9f};


#endif  // TUNING_H_
"""
        )

except IOError:
    print("Cannot open file!")
    exit(1)

except Exception as e:
    print("Error occurred:", e)
    exit(1)
