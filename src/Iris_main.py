import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from NN import Neuron
from GA import ComputeNextGeneration_Storn, ComputeNextGeneration_Smith, ComputeBestKid

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

#define Ackley's function(2D--> 2 variances)
#[(-5<X<5), (-5<Y<5)]
def ComputeFitness_Ackley(DNA):
	X = DNA[0]
	Y = DNA[1]
	#math.exp(1) = 2.71828歐拉常數Euler constant
	fitness = 20 + math.exp(1) - 20 * math.exp(-0.2*(math.sqrt((X**2 + Y**2)/2) - math.exp(0.5 * (math.sin(2 * math.pi * X) + math.cos(2 * math.pi * Y) ))
	return fitness


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
			if (predictions[j] > max_pre_val):
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


#========================================
#主程式開始
#========================================
NO_KIDS = 50

NO_VAR = 35 + 8
#NO_VAR = D(dimension 參數維度)
NO_GEN = 1000

#Hyperparameters setting
MF = 0.7
#mutant factor (0~2)
CR = 0.9
#crossover constant (0~1)

FR = 1.0

SIGMA = 1.0
#========================================
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
kid_dna = -0.5 + 1.0*np.random.rand(NO_KIDS, NO_VAR)
kid_fitness = np.empty(NO_KIDS)

for i in range(NO_KIDS):
	kid_fitness[i] = ComputeFitness(kid_dna[i,:], training_data)
Best_kid = ComputeBestKid(kid_fitness)

#GA/DE algorithm 抓取該世代中fitness值最佳(fitness 最小者，也就是total error 最低)的那一個個體計算並紀錄準確率
#因為輸入ComputeAccuracy function的dna為dna[i,:]，其中i為擁有最佳fitness個體的index
#===================================================
#User input
SOLVER = input("Choose your solver(Smith or Storn): ")
if SOLVER not in ["Smith", "Storn"]:
	print("Invalid input; Smith algorithm will be used by default ")
	SOLVER = "Smith"
#===================================================
history_accuracy = np.empty(NO_GEN)
history_fitness = np.empty(NO_GEN)

for gen in range(NO_GEN):
	if (SOLVER == "Smith"):
		kid_dna, kid_fitness = ComputeNextGeneration_Smith(kid_dna, kid_fitness, Best_kid, FR, SIGMA, ComputeFitness, training_data)
		Best_kid = ComputeBestKid(kid_fitness)
	elif (SOLVER == "Storn"):
		kid_dna, kid_fitness = ComputeNextGeneration_Storn(kid_dna, kid_fitness, MF, CR, ComputeFitness, training_data)
		Best_kid = ComputeBestKid(kid_fitness)
		#Best_kid is the index which has the best fitness


	#calculate accuracy for each generation
	current_accuracy = ComputeAccuracy(kid_dna[Best_kid,:])
	history_accuracy[gen] = current_accuracy
	#history_fitness[gen] = ComputeFitness(kid_dna[Best_kid,:], training_data)
	history_fitness[gen] = kid_fitness[Best_kid]
	#====================================================
	#gen的值從0開始到(NO_GEN-1)結束，應此應該寫(gen + 1)
	print(f"Epoch{gen + 1}: -BestAccuracy : {history_accuracy[gen]*100:.2f}% (Error : {history_fitness[gen]:.4f})")

Final_Accuracy = ComputeAccuracy(kid_dna[Best_kid,:])
Best_Total_Error = kid_fitness[Best_kid]
print("==============final report===============")
print(f"Select SOLVER : {SOLVER}")
print(f"Best Total Error: {Best_Total_Error:.4f}")
print(f"Final Accuracy: {Final_Accuracy*100:.2f}%")

#==========================================
#Smith and Storn
#==========================================

#fig,(ax1,ax2) = plt.subplots(nrows = 2, ncols = 1, figsize = (10,10))
fig,(ax1,ax2) = plt.subplots(2, 1, figsize = (10, 10))
ax1.semilogy(history_fitness, 'r-o', markersize = 3)
#semilogy
#history_fitness 被認定為縱軸（Y 軸）的資料
#y軸為對數刻度
#橫軸（X 軸）則由 Matplotlib 自動以 history_fitness 元素的索引位置來填入
ax1.set(xlabel = 'Generation Number', ylabel = 'Total Error(Fitness)', title = f'[{SOLVER}] GA Training Neural Network - Total Error Descent')
ax1.grid(True)
#開啟網格

ax2.plot(history_accuracy, 'b-')
ax2.set(xlabel = 'Epoch', ylabel = 'Testing Accuracy', title = f'[{SOLVER}] GA Training Neural Network Accuracy on Iris')
ax2.set_ylim(0, 1.05)
ax2.grid(True)
plt.tight_layout()
plt.show()


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
