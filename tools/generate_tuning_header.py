import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, classification_report
import pickle

CSV_FILE = sys.argv[1]
HEADER_FILE = sys.argv[2]
MODEL_FILE = HEADER_FILE.replace('.h', '_poly.pkl')

# Read output from tuning program
df = pd.read_csv(CSV_FILE)

# CSV has columns: M/N, N, Combn, Glynn, Ryser, Fastest
# Fastest column: 0 = Ryser, 1 = Combinatoric, 2 = Glynn
# We only want to use Ryser (0) and Glynn (2) data
# Filter out Combinatoric (1)
tuning = df[df['Fastest'].isin([0, 2])].copy()

# Remap Fastest: 0 (Ryser) stays 0, 2 (Glynn) becomes 1
tuning['Fastest'] = tuning['Fastest'].replace({2: 1})

print(f"Total configurations: {len(tuning)}")
print(f"Algorithm distribution:")
print(tuning['Fastest'].value_counts().sort_index())

features = tuning[["N", "M/N"]].values
labels = tuning["Fastest"].values.astype(int)

# Split into train and test sets
label_counts = pd.Series(labels).value_counts()
can_stratify = all(label_counts >= 2)

if can_stratify:
    x_train, x_test, y_train, y_test = train_test_split(
        features, labels, test_size=0.1, random_state=42, stratify=labels
    )
    print("Using stratified split")
else:
    print(f"Warning: Cannot stratify - some classes have < 2 samples: {label_counts.to_dict()}")
    x_train, x_test, y_train, y_test = train_test_split(
        features, labels, test_size=0.1, random_state=42
    )

# Train polynomial logistic regression
print("\nTraining polynomial logistic regression with grid search...")

pipeline = Pipeline([
    ('poly_features', PolynomialFeatures()),
    ('scaler', StandardScaler()),
    ('logistic', LogisticRegression(max_iter=5000, solver='lbfgs'))
])

param_grid = {
    'poly_features__degree': [1, 2, 3, 4, 5],
    'poly_features__include_bias': [True],
    'logistic__C': [0.01, 0.1, 1, 10, 100]
}

grid_search = GridSearchCV(pipeline, param_grid, cv=5, scoring='accuracy', n_jobs=-1, verbose=1)
grid_search.fit(x_train, y_train)

best_model = grid_search.best_estimator_
best_degree = grid_search.best_params_['poly_features__degree']
best_C = grid_search.best_params_['logistic__C']

print(f"\nBest parameters:")
print(f"  Polynomial degree: {best_degree}")
print(f"  C: {best_C}")
print(f"  Best CV score: {grid_search.best_score_:.4f}")

y_pred = best_model.predict(x_test)
test_accuracy = accuracy_score(y_test, y_pred)
print(f"\nTest accuracy: {test_accuracy:.4f}")

unique_test_classes = np.unique(y_test)
class_names = ['Ryser', 'Glynn']
present_class_names = [class_names[i] for i in unique_test_classes]

print("\nClassification report:")
print(classification_report(y_test, y_pred, labels=unique_test_classes, target_names=present_class_names, zero_division=0))

# Save the model
with open(MODEL_FILE, 'wb') as f:
    pickle.dump(best_model, f)
print(f"\nModel saved to: {MODEL_FILE}")

# Extract model parameters
poly_transformer = best_model.named_steps['poly_features']
scaler = best_model.named_steps['scaler']
log_reg = best_model.named_steps['logistic']

scaler_mean = scaler.mean_
scaler_scale = scaler.scale_
coefficients = log_reg.coef_ 
intercepts = log_reg.intercept_ 

n_features = len(scaler_mean)
n_classes = len(log_reg.classes_) 

# For binary classification, sklearn stores 1 set of coefs (for class 1), but we need both classes
# sklearn: decision_function > 0 means class 1, < 0 means class 0
# We want: argmax of [decision_class0, decision_class1]
# So: decision_class0 = -decision_function, decision_class1 = +decision_function
# But sklearn labels are [0, 1], and we want class 0 = Ryser (0), class 1 = Glynn (1)
# Check class mapping
if n_classes == 2 and coefficients.shape[0] == 1:
    # sklearn stores coefficients for the positive class (higher label)
    # For classes [0, 1], positive class is 1
    # decision_function > 0 → class 1, < 0 → class 0
    # To use argmax: need decision[0] for class 0, decision[1] for class 1
    # decision[0] = -decision_function, decision[1] = +decision_function
    coefficients = np.vstack([-coefficients[0], coefficients[0]])
    intercepts = np.array([-intercepts[0], intercepts[0]])

print(f"\nModel details:")
print(f"  Polynomial degree: {best_degree}")
print(f"  Number of features: {n_features}")
print(f"  Number of classes: {n_classes}")

# Generate C++ header with just the learned parameters
try:
    with open(HEADER_FILE, "w") as file_ptr:
        means_str = ", ".join([f"{m:.15e}" for m in scaler_mean])

        scales_str = ", ".join([f"{s:.15e}" for s in scaler_scale])

        coefs_list = []
        for c in range(n_classes):
            for f in range(n_features):
                coefs_list.append(f"{coefficients[c, f]:.15e}")
        coefs_str = ", ".join(coefs_list)

        intercepts_str = ", ".join([f"{i:.15e}" for i in intercepts])

        file_ptr.write(
            f"""/* Copyright 2024 QC-Devs (GPLv3) */
/* Auto-generated tuning parameters for polynomial logistic regression
 *
 * Model: Polynomial Logistic Regression (degree {best_degree})
 * Train accuracy: {grid_search.best_score_:.4f}
 * Test accuracy: {test_accuracy:.4f}
 *
 * To regenerate: python tools/generate_tuning_header.py <benchmark.csv> <output.h>
 */

#if !defined(permanent_tuning_h_)
#define permanent_tuning_h_

namespace permanent {{

// Model hyperparameters
constexpr int POLY_DEGREE = {best_degree};
constexpr int N_FEATURES = {n_features};
constexpr int N_CLASSES = {n_classes};

// Feature scaling parameters (StandardScaler)
constexpr double SCALER_MEAN[N_FEATURES] = {{ {means_str} }};
constexpr double SCALER_SCALE[N_FEATURES] = {{ {scales_str} }};

// Logistic regression coefficients [N_CLASSES * N_FEATURES]
// Ordered as: [class0_feature0, class0_feature1, ..., class1_feature0, ...]
constexpr double COEFFICIENTS[N_CLASSES * N_FEATURES] = {{ {coefs_str} }};

// Logistic regression intercepts [N_CLASSES]
constexpr double INTERCEPTS[N_CLASSES] = {{ {intercepts_str} }};

}}  // namespace permanent

#endif  // permanent_tuning_h_
"""
        )
    print(f"\nHeader file written to: {HEADER_FILE}")

except IOError:
    print("Cannot open file!")
    exit(1)

except Exception as e:
    print("Error occurred:", e)
    exit(1)
