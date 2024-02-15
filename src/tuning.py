import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.metrics import accuracy_score
import tuning

HEADER_FILE = "permanent/tuning.h"


#load the data
pyTrainDataCombn = tuning.get_pyTrainDataCombn()
pyTestDataCombn = tuning.get_pyTestDataCombn()
pyTrainDataGlynn = tuning.get_pyTrainDataGlynn()
pyTestDataGlynn = tuning.get_pyTestDataGlynn()
pyRyserParam4 = tuning.get_pyRyserParam4()


print(f"Py train data combn: \n{pyTrainDataCombn}")




df = pd.read_csv('fast_permanent.csv', usecols = [' Size', 'M/N', ' Fastest Algorithm'])

# Update label columns for ML
df.rename(columns={' Size': 'x', 'M/N': 'y', ' Fastest Algorithm': 'target'}, inplace=True)

# Make dataset dual class only
update_target = df['target'] == ' Ryser'
df.loc[update_target, 'target'] = ' Combn'

# Update classes to -1/1, Combn = -1
update_combn = df['target'] == ' Combn'
df.loc[update_combn, 'target'] = -1

update_glynn = df['target'] == ' Glynn'
df.loc[update_glynn, 'target'] = 1

# Change feature type to float
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

### Train a linear model - add a hard margin
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
#plt.show()

# Check the accuracy of the model
predictions_linear = linear_model.predict(x_test)
accuracy_linear = accuracy_score(y_test, predictions_linear)
print("Linear Kernel with Hard Margin\nAccuracy (normalized): " + str(accuracy_linear))

# Get the coefficients and bias
coefficients = linear_model.coef_[0]
bias = linear_model.intercept_[0]

# Print the equation of the decision boundary
equation = f"Decision Boundary Equation: {coefficients[0]} * x1 + {coefficients[1]} * x2 + {bias} = 0"
print(equation)

# Write a header file with constants defined as macros 
param_1 = coefficients[0]  # Replace with your actual value
param_2 = coefficients[1]  # Replace with your actual value
param_3 = bias  # Replace with your actual value
param_4 = ryser_limit

try:
    with open(HEADER_FILE, "w") as file_ptr:
        file_ptr.write("#ifndef PERMANENT_TUNING_H\n#define PERMANENT_TUNING_H\n\n\n#define PARAM_1 %.9f\n#define PARAM_2 %.9f\n#define PARAM_3 %.9f\n\n\n#endif /* PERMANENT_TUNING_H */\n" % (param_1, param_2, param_3))
except IOError:
    print("Cannot open file!")
    exit(-1)
except Exception as e:
    print("Error occurred:", e)
    exit(-1)



