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
	rnbr = np.empty(NO_KIDS, dtype = int)
	randb = np.empty(NO_VAR)
	mutant_vector = np.empty(NO_VAR)
	new_fitness = FITNESS.copy()
	new_dna = DNA.copy()
	for i in range(NO_KIDS):
		ParentA = i
		ParentB = i
		ParentC = i
		while((ParentA == ParentB) or (ParentB == ParentC) or (ParentA == ParentC) or (ParentA == i) or (ParentB == i) or (ParentC == i)):
			ParentA = np.random.randint(NO_KIDS)
			ParentB = np.random.randint(NO_KIDS)
			ParentC = np.random.randint(NO_KIDS)
		rnbr[i] = np.random.randint(0, NO_VAR)
		for j in range(NO_VAR):
			randb[j] = np.random.rand()
			if ((randb[j] <= CR) or (j == rnbr[i])):
				mutunt_vector[j] = DNA[ParentA, j] + MF * (DNA[ParentB, j] - DNA[ParentC, j])
			else:
				mutant_vector[j] = DNA[i, j]
		mutant_fitness = ComputeFitness(mutant_vector, traning_data)
		if (mutant_fitness < FITNESS[i]):
			new_fitness[i] = mutant_fitness
			new_dna[i,:] = mutant_vector
	return new_dna, new_fitness


def ComputeNextGeneration_Smith(DNA, FITNESS, FR, SIGMA, ComputeFitness, training_data):
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
		if (trial_fitness < FITNESS[i]):
			new_fitness[i] = trialfitness
			new_dna[i, :] = trial_dna

	return new_dna, new_fitness
