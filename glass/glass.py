import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler



# Decide to use batching or not
USE_BATCHING = True


'''
Load features from file
'''


#feathures = np.delete(data,10,axis=1) 只刪除了最後一欄，即第11欄


data = np.loadtxt('./data/glass.data', delimiter=',')
#Glass Identification Data Set has 214 samples, each with 1ID, 9 features and 1 label (the 11th column)

features = data[:, 1:10].copy() 
#1:10表示只取索引 1開始到索引10之前結束(不包含索引10)，即第二欄到第十欄之間的9欄
#索引0是ID，索引1到9是特徵，：表示取所有列

'''
Load labels from file
'''

#labels = data[:, 10].copy()  #只取第11欄即最後一欄作為標籤
labels = (data[:, 10].copy() - 1).astype(int)
#Glass dataset lables are from 1 to 7(with 4 being empty), we convert them to 0 to 6 for tf.keras

#scaler = StandardScaler()
#features = scaler.fit_transform(features)

# Create a tf.data.Dataset
dataset = tf.data.Dataset.from_tensor_slices((features, labels))

# Iterate and print elements
for element in dataset:
    print(element)

# Split into train / test
num_samples = features.shape[0]
indices = np.arange(num_samples)
np.random.seed(100)
np.random.shuffle(indices)
train_split = int(0.7 * num_samples)
train_idx, test_idx = indices[:train_split], indices[train_split:]

x_train, y_train = features[train_idx], labels[train_idx]
x_test, y_test = features[test_idx], labels[test_idx]



# If we are using batching
if USE_BATCHING:
    batch_size = 16
    train_ds = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(100).batch(batch_size)
    test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(batch_size)



# Model: input -> one hidden dense -> output (7 classes)
#architecture 9-19-7was used in the paper"An ensemble of differential evolution and Adam for training feed-forward neural networks"
model = tf.keras.Sequential([
    tf.keras.layers.Input(shape=(9,)),
    tf.keras.layers.Dense(19, activation='sigmoid'),
    tf.keras.layers.Dense(7, activation='softmax')
])

model.compile(
    #optimizer='adam',
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
    loss=tf.keras.losses.SparseCategoricalCrossentropy(),
    metrics=['accuracy']
)

model.summary()


print("================ Model Construction Complete. Starting training ============")

if USE_BATCHING:
    history = model.fit(train_ds, epochs=1000, validation_data=test_ds)
else:
    history = model.fit(x_train, y_train, epochs=1000, validation_data=(x_test, y_test))

# TODO: Save the model to file to we can load it later without retraining
model.save('glass_model.keras')

print("================Training complete.============")
for i, layer in enumerate(model.layers):
    weights, biases = layer.get_weights()
    print(f"\n===== Layer {i} =====")
    print("Weights:\n", weights)
    print("Biases:\n", biases)

print("Evaluating on test data...")

# Evaluate
if USE_BATCHING:
    loss, acc = model.evaluate(test_ds)
else:
    loss, acc = model.evaluate(x_test, y_test)