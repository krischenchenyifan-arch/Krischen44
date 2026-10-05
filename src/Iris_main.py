import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd

'''
#inputs, biases and weights for testing Run_Network function in class Neuron.
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
'''

'''
#I have moved this to NN.py
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
		#print(Neuron.Run_Network(inputs, weights, biases))
'''

#cost function
def ComputeFitness(DNA, training_data):
	weights = []
	biases = [0, 0, 0, 0]
	for i in range(35):
		weights.append(DNA[i])
	for i in range(8):
		biases.append(DNA[35 + i])
	total_error = 0

	for data in training_data:
		predictions = Neuron.Run_Network(data[0], weights, biases)
		#predictions = Neuron.Run_Network(data['inputs'], weights, biases)
		sample_error = 0
		for i in range(len(predictions)):
			sample_error += (predictions[i] - data[1][i])**2
			#sample_error += (predictions[i] - data['expected'][i])**2
		total_error += sample_error
	return total_error/len(training_data)
	#MSE(Mean Squared Error)
	#neuron 0~3 bias = 0, total 8 biases
	#for i in range(8) i 從0開始記數，所以設定i + 4

 #Accuracy

def ComputeAccuracy(DNA):
	weights = []
	biases = [0, 0, 0, 0]
	for i in range(35):
		weights.append(DNA[i])
	for i in range(8):
		biases.append(DNA[i + 35])
	correct_count = 0
	for data in training_data:
		predictions = Neuron.Run_Network(data[0], weights, biases)
		expected = data[1]
		best_pre_index = 0
		max_pre_val = predictions[0]
		for j in range(1,len(predictions)):
			if (predictons[j] > max_pre_val):
				best_pre_index = j
				max_pre_val = predictions[j]
		best_data_index = 0
		max_data_val = expected[0]
		for k in range(1,len(expected)):
			if (expected[k] > max_data_val):
				best_data_index = k
				max_data_val = expected[k]
		if (best_pre_index == best_data_index):
			correct_count += 1

	return (correct_count/len(training_data))



'''
def ComputeNextGeneration_Storn(DNA, FITNESS, MF, CR, ComputeFitness, training_data):
	NO_KIDS, NO_VAR = DNA.shape
	#若有寫這行就不需要NO_VAR, NO_VAR在函數裡
	mutant_vector = np.empty(NO_VAR)
	randb = np.empty(NO_VAR)
	rnbr = np.empty(NO_KIDS, dtype=int)
	New_dna = DNA.copy()
	New_fitness = FITNESS.copy()
	#gen_val = 0
	#我把initialization放在主程式中
	#現在函數的DNA就是initialization後的初始DNA 
	#while gen_val < NO_GEN:
	#for gen in range(NO_GEN):
	for i in range(NO_KIDS):
		ParentA = i
		ParentB = i
		ParentC = i
		while ((ParentA == ParentB) or (ParentA == ParentC) or (ParentB == ParentC) or (ParentA == i) or (ParentB == i) or (ParentC == i)):
			ParentA = np.random.randint(NO_KIDS)
			ParentB = np.random.randint(NO_KIDS)
			ParentC = np.random.randint(NO_KIDS)
			#np.random.randint(0,44) 與np.random.randint(44)意義相同------->隨機生成0～43之間(含0和43)的整數
			#ParentA, ParentB, ParentC皆為數字(第幾個個體)
		rnbr[i] = np.random.randint(0, NO_VAR)
		#論文中寫到rnbr[i]為random chosen index(1,2,...,D)
		#原本我寫rnbr[i] = np.random.randint(1, NO_VAR + 1)，若rnbr[i]選到43，但j最大只到42，導致j == rnbr[i]不可能發生
		for j in range(NO_VAR):
			#Mutation
			#mutant_vector[j] = DNA[ParentA, j] + MF*(DNA[ParentB, j] - DNA[ParentC, j])
			#Crossover
			randb[j] = np.random.rand()
			#randb[j] is the jth evaluation of a uniform random number generator with outcome between 0 and 1 i.e. [0,1]
			#rnbr[i] = np.random.randint(0, NO_VAR + 1) 
			if ((randb[j] <= CR) or (j == rnbr[i])):
				mutant_vector[j] = DNA[ParentA, j]  + MF*(DNA[ParentB, j] - DNA[ParentC, j])
			else:
				mutant_vector[j] = DNA[i,j] 
		mutant_fitness = ComputeFitness(mutant_vector, training_data)
		#target_fitness[i] = Neuron.ComputeFitness(DNA[i,:] , training_data)
		#target_fitness就是initial_dna(DNA)算出的fitness(函數中輸入的FITNESS)，移動到主程式設定
		if (mutant_fitness < FITNESS[i]):
			New_fitness[i] = mutant_fitness 
			New_dna[i,:] = mutant_vector
	return New_dna, New_fitness 
'''
#========================================
#主程式開始
#========================================
NO_KIDS = 20
#一代有20個個體
NO_VAR = 35 + 8
#NO_VAR = D(dimension 參數維度)
NO_GEN = 400
MF = 1.0
#mutant factor (0~2)
CR = 0.5
#crossover constant (0~1)

df = pd.read_csv('iris_lazy.data', header = None)

training_data = []

for row in df.values.tolist():
	inputs = row[:4]
	#[:4]不包含4即0,1,2,3
	#每一筆資料中前4項為input，最後一項為output
	label = int(row[4])

	if label == 0:
		expected = [1.0, 0.0, 0.0]
	elif label == 1:
		expected = [0.0, 1.0, 0.0]
	else:
		expected = [0.0, 0.0, 1.0]
	#training_data.append({inputs}, {expected})
	#list的append只能一次接受一個參數
	training_data.append([inputs, expected])
	#list of list
#initialization(initial population setting)
initial_dna = -2 + 4*np.random.rand(NO_KIDS, NO_VAR)

'''
#20個(NO_KIDS)初代個體各自的fitness
first_gen_fitness = []
#Kids = [0] * NO_KIDS
for i in range(NO_KIDS):
	#Kid[i] = Neuron.ComputeFitness(initial_dna[i,:], training_data)
	Kid = Neuron.ComputeFitness(initial_dna[i,:], training_data)
	first_gen_fitness.append(Kid)
print(first_gen_fitness)

'''
