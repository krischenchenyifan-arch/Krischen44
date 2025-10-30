import tensorflow as tf
from tensorflow import keras

print("🚀 TensorFlow Hello World is running!")
# 載入資料集
mnist = keras.datasets.fashion_mnist
(train_images, train_labels), (test_images, test_labels) = mnist.load_data()

# 正規化圖片資料（0~255 → 0~1）
train_images = train_images / 255.0
test_images = test_images / 255.0

# 建立模型
model = keras.Sequential([
    keras.layers.Flatten(input_shape=(28, 28)),
    keras.layers.Dense(128, activation='relu'),
    keras.layers.Dense(10)
])

# 編譯模型
model.compile(optimizer='adam',
              loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
              metrics=['accuracy'])

# 訓練模型
model.fit(train_images, train_labels, epochs=1)

# 評估模型
test_loss, test_acc = model.evaluate(test_images, test_labels, verbose=2)
print('\nTest accuracy:', test_acc)
