import tensorflow as tf
import numpy as np
import pandas as pd

# Decide to use batching or not
USE_BATCHING = True


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


# Split into train / test
num_samples = features.shape[0]
indices = np.arange(num_samples)
np.random.seed(42)
np.random.shuffle(indices)
train_split = int(0.8 * num_samples)
train_idx, test_idx = indices[:train_split], indices[train_split:]

x_train, y_train = features[train_idx], labels[train_idx]
x_test, y_test = features[test_idx], labels[test_idx]

# If we are using batching
if USE_BATCHING:
    batch_size = 16
    train_ds = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(100).batch(batch_size)
    test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(batch_size)



print("================ Data splitting and loading complete. Starting Neural Net (Model) construction ============")

# Model: input -> one hidden dense -> output (3 classes)
model = tf.keras.Sequential([
    tf.keras.layers.Input(shape=(4,)),
    tf.keras.layers.Dense(5, activation='relu'),
    tf.keras.layers.Dense(4, activation='relu'),
    tf.keras.layers.Dense(3, activation='softmax')
])

model.compile(
    optimizer='adam',
    loss=tf.keras.losses.SparseCategoricalCrossentropy(),
    metrics=['accuracy']
)

model.summary()

print("================ Model Construction Complete. Starting training ============")

if USE_BATCHING:
    history = model.fit(train_ds, epochs=200, validation_data=test_ds)
else:
    history = model.fit(x_train, y_train, epochs=200, validation_data=(x_test, y_test))

# TODO: Save the model to file to we can load it later without retraining
# model.save('iris_model.h5')

print("================Training complete.============")

print("Evaluating on test data...")

# Evaluate
if USE_BATCHING:
    loss, acc = model.evaluate(test_ds)
else:
    loss, acc = model.evaluate(x_test, y_test)

print(f"Test loss: {loss:.4f}, Test accuracy: {acc:.4f}")