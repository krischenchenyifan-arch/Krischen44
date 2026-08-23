import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd


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
    
    # hidden layer outputs
    hidden_output = []
    for idx, indices in enumerate(hidden_weight_indices):
        neuron_id = idx + 4
        z = sum(Neuron(inputs[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
        hidden_output.append(Neuron.sigmoid(z))
        
    # output layer output
    final_output = []
    for idx, indices in enumerate(output_weight_indices):
        neuron_id = idx + 9
        z = sum(Neuron(hidden_output[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
        final_output.append(Neuron.sigmoid(z))
        
    return final_output


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
        
    one_data = {
        "inputs": inputs,
        "expected": expected
    }
    training_data.append(one_data)


def ComputeFitness(x):    
    weights_dict = {}
    bias_dict = {0: 0, 1: 0, 2: 0, 3: 0}
    
    for i in range(35):
        weights_dict[i] = x[i]
        
    for i in range(8):
        neuron_id = i + 4
        weights_dna_index = 35 + i
        bias_dict[neuron_id] = x[weights_dna_index]
        
    total_error = 0.0
    
    for data in training_data:
        inputs = data["inputs"]
        expected = data["expected"]
        
        predictions = Run_Network(inputs, weights_dict, bias_dict)
        sample_error = 0.0
        for p, e in zip(predictions, expected):
            squared_diff = (p - e) ** 2
            sample_error += squared_diff
            
        total_error += sample_error
        
    return total_error


def ComputeBestKid(FITNESS):
    NO_KIDS = len(FITNESS)
    BestFitness = 10000.0
    BestIndex = -1
    for i in range(NO_KIDS):
        if (np.absolute(FITNESS[i]) < BestFitness):
            BestFitness = FITNESS[i]
            BestIndex = i
    return BestIndex


def MaxMod(X):
    if (X > 0.5):
        return X
    else:
        return (1.0 - X)


def ComputeNextGeneration(DNA, FITNESS, BESTINDEX, FR, SIGMA):
    NO_KIDS, NO_VAR = DNA.shape
    trial_dna = np.empty(NO_VAR)

    new_dna = DNA.copy()
    new_fitness = FITNESS.copy()

    for i in range(NO_KIDS):
        Parent_A = BESTINDEX
        Parent_B = Parent_A
        Parent_C = Parent_A
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


NO_KIDS = 20
NO_VAR = 35 + 8   # 35 weights + 8 biases = 43 
NO_GEN = 15
FR = 1.0
SIGMA = 1.0

kid_dna = -2 + 4 * np.random.rand(NO_KIDS, NO_VAR)
kid_fitness = np.zeros(NO_KIDS)
history_fitness = np.empty(NO_GEN)

for i in range(NO_KIDS):
    kid_fitness[i] = ComputeFitness(kid_dna[i, :])

BestKid = ComputeBestKid(kid_fitness)
print(f"Initial Best Error: {kid_fitness[BestKid]}")

for gen in range(NO_GEN):
    kid_dna, kid_fitness = ComputeNextGeneration(kid_dna, kid_fitness, BestKid, FR, SIGMA)
    BestKid = ComputeBestKid(kid_fitness)
    history_fitness[gen] = kid_fitness[BestKid]
    print(f"Iteration {gen} - Best Total Error = {kid_fitness[BestKid]}")

fig, ax = plt.subplots(figsize=(10, 6))
ax.semilogy(history_fitness, 'r-o')
ax.set(xlabel='Generation Number', ylabel='Total Error (Fitness)', title='GA Training Neural Network on Iris Dataset')
ax.grid(True)
plt.show()

print("====== FINAL REPORT ========")
print(f"Best Total Error: {kid_fitness[BestKid]}")
