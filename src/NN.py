import math

class Neuron:
	def __init__(self, input_val, weight_val):
		self.input = input_val
		self.weight = weight_val
	'''
	def output(self):
		return self.input*self.weight
	'''


	def sigmoid(z):
		if z < -709:
			return 0.0
		elif z > 709:
			return 1.0
		return 1/(1 + math.exp(-z))
	'''
	def sigmoid(z):
                return 1/(1 + math.exp(-z))
	'''
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
			#z = sum(Neuron(inputs[i], weights[w_idx]).output() for i, w_idx in enumerate(indices)) + biases[neuron_id]
			z = sum(inputs[i]*weights[w_idx] for i, w_idx in enumerate(indices)) + biases[neuron_id]
			hidden_output.append(Neuron.sigmoid(z))

		final_output = []
		for idx, indices in enumerate(output_weight_indices):
			neuron_id = idx + 9
			#z = sum(Neuron(hidden_output[i], weights[w_idx]).output() for i, w_idx in enumerate(indices)) + biases[neuron_id]
			z = sum(hidden_output[i]*weights[w_idx] for i, w_idx in enumerate(indices)) + biases[neuron_id]
			final_output.append(Neuron.sigmoid(z))

		return final_output
