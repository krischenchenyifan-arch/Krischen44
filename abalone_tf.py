import tensorflow as tf
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from tensorflow.keras.callbacks import EarlyStopping

# ==============================
# 1. 載入資料與前處理
# ==============================

# 讀取完整資料（包含 Sex）
data = pd.read_csv(
    './abalone/abalone.data',
    header=None,
    names=[
        'Sex', 'Length', 'Diameter', 'Height',
        'WholeWeight', 'ShuckedWeight', 'VisceraWeight', 'ShellWeight', 'Rings'
    ]
)

# One-hot encode 性別
data = pd.get_dummies(data, columns=['Sex'], drop_first=False)

# 分離特徵與標籤
features = data.drop('Rings', axis=1).values
labels = data['Rings'].values

#  印出每一筆資料的 Rings 值
print("\n================ 所有樣本的 Rings 值 ==================")
for i, r in enumerate(labels):
    print(f"Sample {i+1}: Rings = {int(r)}")
print(f"=======================================================\n共 {len(labels)} 筆樣本。")


# 標準化特徵
scaler = StandardScaler()
features = scaler.fit_transform(features)

print(f" Dataset loaded: {features.shape[0]} samples, {features.shape[1]} features")

# ==============================
# 2. 分割訓練與測試集
# ==============================
num_samples = features.shape[0]
indices = np.arange(num_samples)
np.random.seed(42)
np.random.shuffle(indices)

train_split = int(0.8 * num_samples)
train_idx, test_idx = indices[:train_split], indices[train_split:]

x_train, y_train = features[train_idx], labels[train_idx]
x_test, y_test = features[test_idx], labels[test_idx]

# ==============================
# 3. 建立 Dataset（含 batching）
# ==============================
batch_size = 32
train_ds = tf.data.Dataset.from_tensor_slices((x_train, y_train)).shuffle(200).batch(batch_size)
test_ds = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(batch_size)

# ==============================
# 4. 建立最佳化模型
# ==============================
model = tf.keras.Sequential([
    tf.keras.layers.Input(shape=(features.shape[1],)),
    tf.keras.layers.Dense(128, activation='relu', kernel_regularizer=tf.keras.regularizers.l2(0.001)),
    tf.keras.layers.Dropout(0.2),
    tf.keras.layers.Dense(64, activation='relu', kernel_regularizer=tf.keras.regularizers.l2(0.001)),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dense(1, activation='linear')  # 回歸輸出
])

model.compile(
    optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
    loss='mse',
    metrics=['mae']
)

model.summary()

# ==============================
# 5. 訓練模型（含 EarlyStopping）
# ==============================

early_stop = EarlyStopping(monitor='val_loss', patience=30, restore_best_weights=True)
# patience = 
# restore_best_weights =


print("\n================ Start Training ============")
history = model.fit(
    train_ds,
    epochs=100,
    validation_data=test_ds,
    callbacks=[],
    verbose=2
)

# Source - https://stackoverflow.com/a
# Posted by Onno Kampman, modified by community. See post 'Timeline' for change history
# Retrieved 2025-11-14, License - CC BY-SA 3.0

for layer in model.layers:
    print(layer.get_config(), layer.get_weights())

#     callbacks=[early_stop],


# ==============================
# 6. 評估模型
# ==============================
print("\n================ Evaluation ============")
loss, mae = model.evaluate(test_ds)
print(f"Test MSE: {loss:.4f}, MAE: {mae:.4f}")

# 預測與 R² 計算
y_pred = model.predict(x_test).flatten()
r2 = r2_score(y_test, y_pred)
corr = np.corrcoef(y_test, y_pred)[0, 1]

print(f"R² (Coefficient of Determination): {r2:.4f}")
print(f"Correlation between predicted & actual Rings: {corr:.4f}")

# ==============================
# 7. 顯示部分預測結果
# ==============================
print("\n Compare predicted vs actual (前 10 筆):")
for i in range(10):
    print(f"Sample {i+1}: Predicted = {y_pred[i]:.2f}, Actual = {y_test[i]}")

# ==============================
# 8. 可選：儲存模型
# ==============================
model.save('abalone_best_model.h5')









