import math

input = [5.1, 3.8, 1.6, 0.2]

bias = {
	0: 0, 1: 0, 2: 0, 3: 0,
	4: -0.0260747, 5: 0.339481, 6: 0.542879, 7: -0.729046, 8: 0.806341,
	9: 0.4566, 10: 0.475939, 11: -0.564546
}

weight = {
	0: 7.52436, 1: 8.84869, 2: -19.8074, 3: -3.54702, 4: 5.22944,
	5: -2.48485, 6: 6.53303, 7: -14.6117, 8: -22.8034, 9: -4.34393,
	10: 5.99496, 11: 0.207013, 12: -1.11588, 13: -6.96244, 14: -18.3058,
	15: -5.75681, 16: -5.80779, 17: -0.535812, 18: -7.00339, 19: -8.97093,
	20: 12.2787, 21: -7.31267, 22: 1.44651,
	23: -4.53042, 24: 7.34816, 25: -5.36046,
	26: -0.124184, 27: 5.96869, 28: -5.73783,
	29: -2.75999, 30: -12.3828, 31: 0.114309,
	32: 4.90098, 33: -10.149, 34: -8.4269
}

# 隱藏層神經元 (4, 5, 6, 7, 8) 分別對應的權重編號
hidden_weight_indices = [
    [0, 5, 10, 15],  # Neuron 4 用的權重 key
    [1, 6, 11, 16],  # Neuron 5 用的權重 key
    [2, 7, 12, 17],  # Neuron 6 用的權重 key
    [3, 8, 13, 18],  # Neuron 7 用的權重 key
    [4, 9, 14, 19]   # Neuron 8 用的權重 key
]

output_weight_indices = [
    [20, 23, 26, 29, 32],
    [21, 24, 27, 30, 33],
    [22, 25, 28, 31, 34]
]


class Neuron:
	def __init__(self, input, weight):
		self.input = input
		self.weight = weight
	def output(self):
		z = self.input*self.weight
		return z
	@staticmethod
	def sigmoid(a):
		return 1/(1 + math.exp(-a))

#Calculate hidden layer neurons' outputs
hidden_output = []
for idx, indices in enumerate(hidden_weight_indices):
	neuron_id = idx + 4
	z = sum(Neuron(input[i], weight[w_idx]).output() for i, w_idx in enumerate(indices)) + bias[neuron_id]
	final_value = Neuron.sigmoid(z)
	hidden_output.append(final_value)
	print(f"Neuron{neuron_id} output is {final_value}")	
print(hidden_output)
print(f"======start calculate output layer outputs======")

final_output = []
for idx, indices in enumerate(output_weight_indices):
	neuron_id = idx + 9
	z = sum(Neuron(hidden_output[i], weight[w_idx]).output() for i,w_idx in enumerate(indices)) + bias[neuron_id]
	final_value = Neuron.sigmoid(z)
	final_output.append(final_value)
	print(f"Neuron{neuron_id} output is {final_value}")
print(final_output)
