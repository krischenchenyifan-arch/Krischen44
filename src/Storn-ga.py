import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

inputs = [5.1, 3.8, 1.6, 0.2]

biases = [
	 0, 0, 0, 0,
	 -0.0260747, 0.339481, 0.542879, -0.729046, 0.806341,
	 0.4566, 0.475939, -0.564546
]

weights = [
	 7.52436, 8.84869, -19.8074, -3.54702, 5.22944,
	 -2.48485, 6.53303, -14.6117, -22.8034, -4.34393,
	 5.99496, 0.207013, -1.11588, -6.96244, -18.3058,
	 -5.75681, -5.80779, -0.535812, -7.00339, -8.97093,
	 12.2787, -7.31267, 1.44651,
	 -4.53042, 7.34816, -5.36046,
	 -0.124184, 5.96869, -5.73783,
	 -2.75999, -12.3828, 0.114309,
	 4.90098, -10.149, -8.4269
]


class Neuron:
	def __init__(self, input_val, weight_val):
		self.input = input_val
		self.weight = weight_val

	def output(self):
		return self.input*self.weight

	def sigmoid(z):
		if z < -709:
			return 0.0
		elif z > 709:
			return 1.0
		return 1/(1 + math.exp(-z))

	def Run_Network(inputs, weights, biases):
		hidden_weight_indices = [
	  	[0, 5, 10, 15],
	  	[1, 6, 11, 16],
	  	[2, 7, 12, 17],
	  	[3, 8, 13, 18],
	  	[4, 9, 14, 19]
		] 

		output_weight_indices = [
	  	[20, 23, 26, 29, 32],
	  	[21, 24, 27, 30, 33],
	  	[22, 25, 28, 31, 34]
		]

		hidden_output = []
		for idx, indices in enumerate(hidden_weight_indices):
			neuron_id = idx + 4
			z = sum(Neuron(inputs[i], weights[w_idx]).output() for i, w_idx in enumerate(indices)) + biases[neuron_id]
			hidden_output.append(Neuron.sigmoid(z))

		final_output = []
		for idx, indices in enumerate(output_weight_indices):
			neuron_id = idx + 9
			z = sum(Neuron(hidden_output[i], weights[w_idx]).output() for i, w_idx in enumerate(indices)) + biases[neuron_id]
			final_output.append(Neuron.sigmoid(z))

		return final_output
print(Neuron.Run_Network(inputs, weights, biases))
