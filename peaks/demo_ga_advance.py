import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

# --- 1. 定義神經網路類別與前向傳播函數 ---
class Neuron:
    def __init__(self, input_val, weight_val):
        self.input = input_val
        self.weight = weight_val
    def output(self):
        return self.input * self.weight
    @staticmethod
    def sigmoid(a):
        return 1 / (1 + math.exp(-a))

def Run_Network(inputs, weights_dict, bias_dict):
    hidden_weight_indices = [
        [0, 5, 10, 15],  # Neuron 4
        [1, 6, 11, 16],  # Neuron 5
        [2, 7, 12, 17],  # Neuron 6
        [3, 8, 13, 18],  # Neuron 7
        [4, 9, 14, 19]   # Neuron 8
    ]
    output_weight_indices = [
        [20, 23, 26, 29, 32],
        [21, 24, 27, 30, 33],
        [22, 25, 28, 31, 34]
    ]
    
    # 計算隱藏層輸出 (Neuron 4~8)
    hidden_output = []
    for idx, indices in enumerate(hidden_weight_indices):
        neuron_id = idx + 4
        z = sum(Neuron(inputs[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
        hidden_output.append(Neuron.sigmoid(z))
        
    # 計算輸出層輸出 (Neuron 9~11)
    final_output = []
    for idx, indices in enumerate(output_weight_indices):
        neuron_id = idx + 9
        z = sum(Neuron(hidden_output[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
        final_output.append(Neuron.sigmoid(z))
        
    return final_output


# --- 2. 直接讀取 Iris 資料集檔案 ---
df = pd.read_csv('iris_lazy.data', header=None)

training_data = []
for index, row in df.iterrows():
    inputs = [row[0], row[1], row[2], row[3]]
    label = int(row[4])
    
    if label == 0:
        expected = [1.0, 0.0, 0.0]
    elif label == 1:
        expected = [0.0, 1.0, 0.0]
    else:
        expected = [0.0, 0.0, 1.0]
        
    training_data.append({"inputs": inputs, "expected": expected})


# --- 3. 定義 ComputeFitness 函數 ---
def ComputeFitness(x):
    weights_dict = {}
    bias_dict = {0: 0, 1: 0, 2: 0, 3: 0}
    
    for i in range(35):
        weights_dict[i] = x[i]
    for i in range(8):
        bias_dict[i + 4] = x[35 + i]
        
    total_error = 0.0
    for data in training_data:
        predictions = Run_Network(data["inputs"], weights_dict, bias_dict)
        sample_error = sum((p - e) ** 2 for p, e in zip(predictions, data["expected"]))
        total_error += sample_error
        
    return total_error


# --- 4. 輔助函數：計算 Accuracy ---
def ComputeAccuracy(dna):
    weights_dict = {i: dna[i] for i in range(35)}
    bias_dict = {0:0, 1:0, 2:0, 3:0}
    for i in range(8):
        bias_dict[i + 4] = dna[35 + i]
        
    correct_count = 0
    for data in training_data:
        predictions = Run_Network(data["inputs"], weights_dict, bias_dict)
        # 用純 Python 的 max 找出最大機率的索引
        best_pred_idx, _ = max(enumerate(predictions), key=lambda x: x[1])
        best_true_idx, _ = max(enumerate(data["expected"]), key=lambda x: x[1])
        
        if best_pred_idx == best_true_idx:
            correct_count += 1
            
    return correct_count / len(training_data)


# --- 5. 基因演算法核心函數 ---
def ComputeBestKid(FITNESS):
    BestFitness = 10000.0
    BestIndex = -1
    for i in range(len(FITNESS)):
        if (np.absolute(FITNESS[i]) < BestFitness):
            BestFitness = FITNESS[i]
            BestIndex = i
    return BestIndex

def MaxMod(X):
    return X if X > 0.5 else (1.0 - X)

def ComputeNextGeneration(DNA, FITNESS, BESTINDEX, FR, SIGMA):
    NO_KIDS, NO_VAR = DNA.shape
    trial_dna = np.empty(NO_VAR)
    new_dna = DNA.copy()
    new_fitness = FITNESS.copy()

    for i in range(NO_KIDS):
        Parent_A = BESTINDEX
        Parent_B, Parent_C = Parent_A, Parent_A
        while ((Parent_A == Parent_B) or (Parent_A == Parent_C) or (Parent_B == Parent_C)):
            Parent_B = np.random.randint(NO_KIDS)
            Parent_C = np.random.randint(NO_KIDS)

        for j in range(NO_VAR):
            Rf = FR * MaxMod(np.random.rand())
            Rnf = SIGMA * np.random.randn()
            trial_dna[j] = Rf * DNA[Parent_A, j] + (1.0 - Rf) * DNA[Parent_B, j] + Rnf * (DNA[Parent_B, j] - DNA[Parent_C, j])

        trial_fitness = ComputeFitness(trial_dna)
        if (trial_fitness < FITNESS[i]):
            new_fitness[i] = trial_fitness
            new_dna[i, :] = trial_dna

    return new_dna, new_fitness


# --- 6. 主程式執行序列 ---
NO_KIDS = 20
NO_VAR = 35 + 8   
NO_GEN = 15      # 為了讓曲線有明顯爬升，可將世代數調高（例如 100 或 400）
FR = 1.0
SIGMA = 1.0

kid_dna = -2 + 4 * np.random.rand(NO_KIDS, NO_VAR)
kid_fitness = np.zeros(NO_KIDS)
history_accuracy = np.empty(NO_GEN) # 改用矩陣記錄每一代的 Accuracy

for i in range(NO_KIDS):
    kid_fitness[i] = ComputeFitness(kid_dna[i, :])

BestKid = ComputeBestKid(kid_fitness)

# 世代演化
for gen in range(NO_GEN):
    kid_dna, kid_fitness = ComputeNextGeneration(kid_dna, kid_fitness, BestKid, FR, SIGMA)
    BestKid = ComputeBestKid(kid_fitness)
    
    # 計算並記錄當代最佳小孩的 Accuracy (0 ~ 1 之間)
    current_best_accuracy = ComputeAccuracy(kid_dna[BestKid, :])
    history_accuracy[gen] = current_best_accuracy
    
    print(f"Generation {gen} - Best Accuracy = {current_best_accuracy * 100:.2f}%")

# --- 7. 繪製類似截圖的 Accuracy 趨勢圖 ---
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(history_accuracy, 'b-') # 用藍色實線繪製
ax.set(xlabel='Epoch', ylabel='Testing Accuracy', title='GA Training Neural Network Accuracy')
ax.set_ylim(0, 1.05) # 設定 Y 軸從 0 到 1.05
ax.grid(True)
plt.show()

print("====== FINAL REPORT ========")
print(f"Final Best Accuracy: {history_accuracy[-1] * 100:.2f}%")
