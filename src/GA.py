import numpy as np
import copy

def MaxMod(X):
	return max(X, 1 - X)
'''
def MaxMod(X):
	if (X > 0.5):
		return X
	else:
		return (1.0 - X)
'''
def ComputeBestKid(FITNESS):
	NO_KIDS = len(FITNESS)
	BestFitness = 10000.0
	BestIndex = -1
	for i in range(NO_KIDS):
		if (np.abs(FITNESS[i]) < BestFitness):
			BestFitness = np.abs(FITNESS[i])
			BestIndex = i
	return BestIndex
def ComputeNextGeneration_Storn(DNA, FITNESS, MF, CR, ComputeFitness):
#def ComputeNextGeneration_Storn(DNA, FITNESS, MF, CR, ComputeFitness, training_data):
	NO_KIDS, NO_VAR = DNA.shape
	#若有寫這行就不需要NO_VAR, NO_VAR在函數裡
	rnbr = np.empty(NO_KIDS, dtype = int)
	randb = np.empty(NO_VAR)
	mutant_vector = np.empty(NO_VAR)
	#new_fitness = copy.deepcopy(FITNESS)
	#new_dna = copy.deepcopy(DNA)
	new_fitness = FITNESS.copy()
	new_dna = DNA.copy()
	#new_fitness = FITNESS
	#new_dna = DNA
	#gen_val = 0
	#我把initialization放在主程式中
	#現在函數的DNA就是initialization後的初始DNA
	#while gen_val < NO_GEN:
	#for gen in range(NO_GEN):
	for i in range(NO_KIDS):
		ParentA = i
		ParentB = i
		ParentC = i
		while((ParentA == ParentB) or (ParentB == ParentC) or (ParentA == ParentC) or (ParentA == i) or (ParentB == i) or (ParentC == i)):
			ParentA = np.random.randint(NO_KIDS)
			ParentB = np.random.randint(NO_KIDS)
			ParentC = np.random.randint(NO_KIDS)
			#np.random.randint(0,44) 與np.random.randint(44)意義相同------->隨機生成0～43之間(含0和43)的整數
			#ParentA, ParentB, ParentC皆為數字(第幾個個體)
		#rnbr[i] = np.random.randint(0, NO_VAR)
		rnbr[i] = np.random.randint(NO_VAR)
		#因為下方for j in range(NO_VAR): j 從0到NO_VAR-1，故rnbr[i] = np.random.randint(0, NO_VAR)，才能在j == rnbr[i]一一對應
		#rnbr[i] = np.random.randint(1, NO_VAR + 1)
		#論文中寫到rnbr[i]為random chosen index(1,2,...,D)
		#原本我寫rnbr[i] = np.random.randint(1, NO_VAR + 1)，若rnbr[i]選到43，但j最大只到42，導致j == rnbr[i]不可能發生
		for j in range(NO_VAR):
			#Mutation
			#Crossover
			randb[j] = np.random.rand()
			#uniform distribution
			#randb[j] is the jth evaluation of a uniform random number generator with outcome between 0 and 1 i.e. [0,1]
			#j = 1,2,...,D
			#對於D維度中的每一個維度j，系統都會抽一個介於[0, 1]之間的亂數randb(j)
			#rnbr[i] = np.random.randint(0, NO_VAR + 1)
			if ((randb[j] <= CR) or (j == rnbr[i])):
				mutant_vector[j] = DNA[ParentA, j] + MF * (DNA[ParentB, j] - DNA[ParentC, j])
			else:
				mutant_vector[j] = DNA[i, j]
		#mutant_fitness = ComputeFitness(mutant_vector, training_data)
		#target_fitness[i] = Neuron.ComputeFitness(DNA[i,:] , training_data)
		#target_fitness就是initial_dna(DNA)算出的fitness(函數中輸入的FITNESS)，移動到主程式設定
		mutant_fitness = ComputeFitness(mutant_vector)
		if (mutant_fitness < FITNESS[i]):
			new_fitness[i] = mutant_fitness
			new_dna[i,:] = mutant_vector
	return new_dna, new_fitness

def ComputeNextGeneration_Smith(DNA, FITNESS, BESTINDEX, FR, SIGMA, ComputeFitness):
#def ComputeNextGeneration_Smith(DNA, FITNESS, BESTINDEX, FR, SIGMA, ComputeFitness, training_data):
	NO_KIDS, NO_VAR = DNA.shape
	trial_dna = np.empty((NO_KIDS, NO_VAR))
	#trial_dna = np.empty(NO_VAR)
	trial_fitness = np.empty(NO_KIDS)
	#new_fitness = copy.deepcopy(FITNESS)
	#new_dna = copy.deepcopy(DNA)
	#new_dna = DNA
	#new_fitness = FITNESS
	new_fitness = FITNESS.copy()
	new_dna = DNA.copy()
	#new_fitness = FITNESS
	#new_dna = DNA
	for i in range(NO_KIDS):
		Parent_A = BESTINDEX
		Parent_B = Parent_A
		Parent_C = Parent_A
		while ((Parent_A == Parent_B) or (Parent_A == Parent_C) or (Parent_B == Parent_C)):
			Parent_B = np.random.randint(NO_KIDS)
			Parent_C = np.random.randint(NO_KIDS)
		for j in range(NO_VAR):
			Gf = FR * MaxMod(np.random.rand())
			Rnf = SIGMA * np.random.randn()
			trial_dna[i, j] = Gf * DNA[Parent_A, j] + (1.0 - Gf) * DNA[Parent_B, j] + Rnf * (DNA[Parent_B, j] - DNA[Parent_C, j])
		#trial_fitness = ComputeFitness(trial_dna, training_data)
		trial_fitness[i] = ComputeFitness(trial_dna[i, :])
		#trial_fitness = ComputeFitness(trial_dna)
		#if (trial_fitness < FITNESS[i]):
		#if (trial_fitness[i] < FITNESS[i]):
		if (np.abs(trial_fitness[i]) < np.abs(FITNESS[i])):
			new_fitness[i] = trial_fitness[i]
			#new_fitness[i] = trial_fitness
			#new_dna[i, :] = trial_dna
			new_dna[i, :] = trial_dna[i, :]
	return new_dna, new_fitness
