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
