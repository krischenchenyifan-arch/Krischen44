import tensorflow as tf
import numpy as np
import pandas as pd


# IRIS data
# IRIS has 4 features
# sepal length	Feature	Continuous		cm	no
# sepal width	Feature	Continuous		cm	no
# petal length	Feature	Continuous		cm	no
# petal width	Feature	Continuous		cm	no
# There are 150 sets
# Each with 4 features



'''
Load features from file
'''

# Example of 3 data sets with 2 inputs per data set
# features = np.array([[1, 2], [3, 4], [5, 6]])

# Try dumb loading the data
# This doesn't work because the iris.data file contains mixed numbers as well as text.
# data = np.loadtxt('./data/iris.data')


# This will work because we got lazy and renamed the last column
data = np.loadtxt('./data/iris_lazy.data', delimiter=',')

# Now we would like to copy this data into a new variable (features) and then delete the last column
# If this copy doesn't work, it might be because python doesn't like performing deep copies of np arrays so lightly

# This is ready now
features = np.delete(data, 4, axis=1)

'''
Load labels from file
'''

# Example of labels - 3 data sets each has a label
# labels = np.array([0, 1, 0])


# Try copying the 4th column from data as our label
labels = data[:, 4].copy()

# print(f"Our labels: {labels}")

# Create a tf.data.Dataset
dataset = tf.data.Dataset.from_tensor_slices((features, labels))

# Iterate and print elements
for element in dataset:
    print(element)