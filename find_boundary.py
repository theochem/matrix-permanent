import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import accuracy_score
from sympy import symbols, simplify


#load the data
df = pd.read_csv('fast_permanent.csv', usecols = [' Size', 'M/N', ' Fastest Algorithm'])

# Update label columns for ML
df.rename(columns={' Size': 'x', 'M/N': 'y', ' Fastest Algorithm': 'target'}, inplace=True)

# Make dataset dual class only
update = df['target'] == ' Ryser'

df.loc[update, 'target'] = ' Combn'

# Update classes to -1/1, Combn = -1
update = df['target'] == ' Combn'
df.loc[update, 'target'] = -1

update = df['target'] == ' Glynn'
df.loc[update, 'target'] = 1

df['target'] = df['target'].astype(float)
features = df[['x', 'y']]
label = df['target']

value_counts = df['target'].value_counts()

# Create train/test split 
size = df.shape[0]
test_size = int(np.round(size * 0.1, 0))

# Split dataset into training and testing sets
x_train = features[:-test_size].values
y_train = label[:-test_size].values
x_test = features[-test_size:].values
y_test = label[-test_size:].values

# Plotting the training set
fig, ax = plt.subplots(figsize=(12, 7))

# Setting the boarders
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adding major gridlines
ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
ax.scatter(features[:-test_size]['x'], features[:-test_size]['y'], color="#8C7298")
plt.show()

# Separate two classes using 2nd-degree polynomial curve
model_poly_hard_margin = svm.SVC(kernel='poly', degree=2, C=10.0)
model_poly_hard_margin.fit(x_train, y_train)

fig, ax = plt.subplots(figsize=(12, 7))

# Setting the boarders

ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

# Create grid to evaluate model
xx = np.linspace(-1, max(features['x']) + 1, len(x_train))
yy = np.linspace(0, max(features['y']) + 1, len(y_train))
YY, XX = np.meshgrid(yy, xx)
xy = np.vstack([XX.ravel(), YY.ravel()]).T
train_size = len(features[:-test_size]['x'])

# Assigning different colors to the classes
colors = y_train
colors = np.where(colors == 1, '#8C7298', '#4786D1')

# Plot the dataset
ax.scatter(features[:-test_size]['x'], features[:-test_size]['y'], c=colors)

# Get the separating hyperplane
Z = model_poly_hard_margin.decision_function(xy).reshape(XX.shape)

# Draw the decision boundary and margins
ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5, linestyles=['--', '-', '--'])

# Highlight support vectors with a circle around them
ax.scatter(model_poly_hard_margin.support_vectors_[:, 0], model_poly_hard_margin.support_vectors_[:, 1], s=100, linewidth=1, facecolors='none', edgecolors='k')
plt.show()

predictions_poly = model_poly_hard_margin.predict(x_test)
accuracy_poly = accuracy_score(y_test, predictions_poly)
print("2nd degree polynomial Kernel\nAccuracy (normalized): " + str(accuracy_poly))

# Get the equation of the decision boundary
support_vectors = model_poly_hard_margin.support_vectors_
alphas = model_poly_hard_margin.dual_coef_.ravel()
intercept = model_poly_hard_margin.intercept_[0]

# Define symbolic variables
x1, x2 = symbols('x1 x2')
degree = 2

# Polynomial kernel parameters
degree = model_poly_hard_margin.degree
coef0 = model_poly_hard_margin.coef0

# Compute the symbolic polynomial kernel
polynomial_kernel = (np.dot(support_vectors, [x1, x2]) + coef0) ** degree

# Combine the symbolic kernel with support vector coefficients
decision_function = sum(alphas[i] * model_poly_hard_margin.support_vectors_[i] * polynomial_kernel[i] for i in range(len(alphas)))

# Add the intercept
decision_function += intercept

# Simplify the expression
simplified_decision_function = simplify(decision_function)

# Print the simplified decision function
print("Simplified Decision Function:")
print(simplified_decision_function)

# def predict_class(x1, x2):
#     # Coefficients obtained from the SVM model - focus on Glynn
#     combn_val = (-17926.1057957553*(x1**2) + 41433.0615998245*(x1*x2) + 1034.96055282324*(x2**2) - 0.688319770680141)
#     glynn_val = 20716.5307999122*(x1**2) + 2069.92110564648*(x1*x2) + 46.5474266407675*(x2**2) - 0.688319770680141

#     print(combn_val)
#     print(glynn_val)
#     return combn_val < glynn_val

# # Test the prediction function
# x1_test, x2_test = 0.5, 12
# predicted_class_test = predict_class(x1_test, x2_test)


# if predicted_class_test is True:
# 	print("Predicted Glynn in the poly case \n")
# else:
# 	print("Predicted Combn in the poly case \n")

# Test a more simple approach - add a hard margin
linear_model = svm.SVC(kernel='linear', C=100.0)
linear_model.fit(x_train, y_train)

fig, ax = plt.subplots(figsize=(12, 7))

# Setting the boarders
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

# Create grid to evaluate model
xx = np.linspace(-1, max(features['x']) + 1, len(x_train))
yy = np.linspace(0, max(features['y']) + 1, len(y_train))
YY, XX = np.meshgrid(yy, xx)
xy = np.vstack([XX.ravel(), YY.ravel()]).T
train_size = len(features[:-test_size]['x'])

# Assigning different colors to the classes
colors = y_train
colors = np.where(colors == 1, '#8C7298', '#4786D1')

# Plot the dataset
ax.scatter(features[:-test_size]['x'], features[:-test_size]['y'], c=colors)

# Get the separating hyperplane
Z = linear_model.decision_function(xy).reshape(XX.shape)

# Draw the decision boundary and margins
ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5, linestyles=['--', '-', '--'])

# Highlight support vectors with a circle around them
ax.scatter(linear_model.support_vectors_[:, 0], linear_model.support_vectors_[:, 1], s=100, linewidth=1, facecolors='none', edgecolors='k')
plt.show()

# Check the accuracy of the model
predictions_linear = linear_model.predict(x_test)
accuracy_linear = accuracy_score(y_test, predictions_linear)
print("Linear Kernel\nAccuracy (normalized): " + str(accuracy_linear))

# Get the coefficients and bias
coefficients = linear_model.coef_[0]
bias = linear_model.intercept_[0]

# Print the equation of the decision boundary
equation = f"Decision Boundary Equation: {coefficients[0]} * x1 + {coefficients[1]} * x2 + {bias} = 0"
print(equation)






