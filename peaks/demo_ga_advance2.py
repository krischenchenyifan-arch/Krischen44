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
        #OverflowError
        if a < -709:
            return 0.0
        elif a > 709:
            return 1.0
        return 1 / (1 + math.exp(-a))

    @staticmethod
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


        hidden_output = []
        for idx, indices in enumerate(hidden_weight_indices):
            neuron_id = idx + 4
            z = sum(Neuron(inputs[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
            hidden_output.append(Neuron.sigmoid(z))

        final_output = []
        for idx, indices in enumerate(output_weight_indices):
            neuron_id = idx + 9
            z = sum(Neuron(hidden_output[i], weights_dict[w_idx]).output() for i, w_idx in enumerate(indices)) + bias_dict[neuron_id]
            final_output.append(Neuron.sigmoid(z))

        return final_output

    @staticmethod
    def ComputeFitness(x, training_data):
	#weights = x[:35]
	#biases = [0.0, 0.0, 0.0, 0.0] + list(x[35:])
        weights_dict = {}
        bias_dict = {0: 0, 1: 0, 2: 0, 3: 0}

        for i in range(35):
            weights_dict[i] = x[i]

        for i in range(8):
            neuron_id = i + 4
            bias_dict[neuron_id] = x[35 + i]

        total_error = 0.0

        for data in training_data:
            predictions = Neuron.Run_Network(data["inputs"], weights_dict, bias_dict)
            sample_error = sum((p - e) ** 2 for p, e in zip(predictions, data["expected"]))
            total_error += sample_error

        return total_error


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


def ComputeAccuracy(dna):
    #weights = dna[:35]
    #biases = [0.0, 0.0, 0.0, 0.0] + list(dna[35:])
    weights_dict = {}
    bias_dict = {0: 0, 1: 0, 2: 0, 3: 0}
    for i in range(35):
       weights_dict[i] = dna[i]
    for i in range(8):
       bias_dict[i + 4] = dna[35 + i]

    correct_count = 0
    for data in training_data:
        predictions = Neuron.Run_Network(data["inputs"], weights_dict, bias_dict)
        expected = data["expected"]

        best_pred_idx = 0
        max_pred_val = predictions[0]
        for i in range(1, len(predictions)):
            if predictions[i] > max_pred_val:
                max_pred_val = predictions[i]
                best_pred_idx = i

        best_true_idx = 0
        max_true_val = expected[0]
        for i in range(1, len(expected)):
            if expected[i] > max_true_val:
                max_true_val = expected[i]
                best_true_idx = i

        if best_pred_idx == best_true_idx:
            correct_count += 1

    return correct_count / len(training_data)


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
    return X if X > 0.5 else (1.0 - X)


def ComputeNextGeneration(DNA, FITNESS, BESTINDEX, FR, SIGMA, training_data):
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

        trial_fitness = Neuron.ComputeFitness(trial_dna, training_data)

        if (trial_fitness < FITNESS[i]):
            new_fitness[i] = trial_fitness
            new_dna[i, :] = trial_dna

    return new_dna, new_fitness


NO_KIDS = 20
NO_VAR = 35 + 8
NO_GEN = 400
FR = 1.0
SIGMA = 1.0

kid_dna = -2 + 4 * np.random.rand(NO_KIDS, NO_VAR)
kid_fitness = np.zeros(NO_KIDS)
history_accuracy = np.empty(NO_GEN)
history_fitness = np.empty(NO_GEN)

for i in range(NO_KIDS):
    kid_fitness[i] = Neuron.ComputeFitness(kid_dna[i, :], training_data)

BestKid = ComputeBestKid(kid_fitness)
print(f"Initial Best Error: {kid_fitness[BestKid]}")

for gen in range(NO_GEN):
    kid_dna, kid_fitness = ComputeNextGeneration(kid_dna, kid_fitness, BestKid, FR, SIGMA, training_data)
    BestKid = ComputeBestKid(kid_fitness)

    current_accuracy = ComputeAccuracy(kid_dna[BestKid, :])
    history_accuracy[gen] = current_accuracy
    history_fitness[gen] = kid_fitness[BestKid]

    print(f"Epoch {gen} - Best Accuracy = {current_accuracy * 100:.2f}% (Error: {kid_fitness[BestKid]:.4f})")

#繪製圖表
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

#Total Error
ax1.semilogy(history_fitness, 'r-o', markersize=3)
ax1.set(xlabel='Generation Number', ylabel='Total Error (Fitness)', title='GA Training Neural Network - Total Error Descent')
ax1.grid(True)

#Testing Accuracy
ax2.plot(history_accuracy, 'b-')
ax2.set(xlabel='Epoch', ylabel='Testing Accuracy', title='GA Training Neural Network Accuracy on Iris')
ax2.set_ylim(0, 1.05)
ax2.grid(True)
plt.tight_layout()
plt.show()

final_accuracy = ComputeAccuracy(kid_dna[BestKid, :])
print("====== FINAL REPORT ========")
print(f"Best Total Error: {kid_fitness[BestKid]:.4f}")
print(f"Final Accuracy: {final_accuracy * 100:.2f}%")
