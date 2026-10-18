import numpy as np


def MaxMod(X):
	if (X > 0.5):
		return X
	else:
		return (1.0 - X)

def ComputeBestKid(FITNESS):
	NO_KIDS = len(FITNESS)
	BestFitness = 10000
	BestIndex = -1
	for i in range(NO_KIDS):
		if (np.absolute(FITNESS[i]) < BestFitness):
			BestFitness = FITNESS[i]
			BestIndex = i
	return BestIndex

def ComputeNextGeneration_Storn(DNA, FITNESS, MF, CR, ComputeFitness, training_data):
	NO_KIDS, NO_VAR = DNA.shape
	#若有寫這行就不需要NO_VAR, NO_VAR在函數裡
	rnbr = np.empty(NO_KIDS, dtype = int)
	randb = np.empty(NO_VAR)
	mutant_vector = np.empty(NO_VAR)
	new_fitness = FITNESS.copy()
	new_dna = DNA.copy()
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
		rnbr[i] = np.random.randint(0, NO_VAR)
		#論文中寫到rnbr[i]為random chosen index(1,2,...,D)
		#原本我寫rnbr[i] = np.random.randint(1, NO_VAR + 1)，若rnbr[i]選到43，但j最大只到42，導致j == rnbr[i]不可能發生
		for j in range(NO_VAR):
			#Mutation
			#Crossover
			randb[j] = np.random.rand()
			#randb[j] is the jth evaluation of a uniform random number generator with outcome between 0 and 1 i.e. [0,1]
			#rnbr[i] = np.random.randint(0, NO_VAR + 1)
			if ((randb[j] <= CR) or (j == rnbr[i])):
				mutant_vector[j] = DNA[ParentA, j] + MF * (DNA[ParentB, j] - DNA[ParentC, j])
			else:
				mutant_vector[j] = DNA[i, j]
		mutant_fitness = ComputeFitness(mutant_vector, training_data)
		#target_fitness[i] = Neuron.ComputeFitness(DNA[i,:] , training_data)
		#target_fitness就是initial_dna(DNA)算出的fitness(函數中輸入的FITNESS)，移動到主程式設定
		#mutant_fitness = ComputeFitness(mutant_vector)
		if (mutant_fitness < FITNESS[i]):
			new_fitness[i] = mutant_fitness
			new_dna[i,:] = mutant_vector
	return new_dna, new_fitness


def ComputeNextGeneration_Smith(DNA, FITNESS, BESTINDEX, FR, SIGMA, ComputeFitness, training_data):
	NO_KIDS, NO_VAR = DNA.shape
	trial_dna = np.empty(NO_VAR)
	new_fitness = FITNESS.copy()
	new_dna = DNA.copy()
	for i in range(NO_KIDS):
		ParentA = BESTINDEX
		ParentB = ParentA
		ParentC = ParentA
		while ((ParentA == ParentB) or (ParentB == ParentC) or (ParentA == ParentC)):
			ParentB = np.random.randint(NO_KIDS)
			ParentC = np.random.randint(NO_KIDS)
		for j in range(NO_VAR):
			Rf = FR*MaxMod(np.random.rand()) 
			Rnf = SIGMA*np.random.randn()
			trial_dna[j] = Rf*DNA[ParentA, j] + (1.0 - Rf)*DNA[ParentB, j] + Rnf*(DNA[ParentB, j] - DNA[ParentC, j])

		trial_fitness = ComputeFitness(trial_dna, training_data)
		#trial_fitness = ComputeFitness(trial_dna)
		if (trial_fitness < FITNESS[i]):
			new_fitness[i] = trial_fitness
			new_dna[i, :] = trial_dna

	return new_dna, new_fitness
