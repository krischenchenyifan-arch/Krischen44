import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from GA import ComputeNextGeneration_Storn, ComputeNextGeneration_Smith


#define 2D Ackley's function(2 variances)
#[(-5 < X < 5), (-5 < Y < 5)]
def Computefitness_Ackley(DNA):
	X = DNA[0]
	Y = DNA[1]
	term1 = -0.2 * math.sqrt((X**2 + Y**2)/2)
	term2 = 0.5 * (math.cos(2*math.pi*X) + math.cos(2*math.pi*Y))
	fitness = 20 + math.exp(1) - 20 * math.exp(term1) - math.exp(term2)
	#math.exp(1) = 2.71828歐拉常數Euler constant
	return fitness

#====================================
#主程式
#====================================
#參數設定
NO_GEN = 50
NO_VAR = 2
NO_KIDS = 20
NUM_RUNS = 20
#The mean best fitness is shown over 1000 unique simulations with each using 
#20 children per generation over 50 generation.
#====================================

for run in range(NUM_RUNS):
	kid_dna = -5 + 5 * np.random.rand(NO_KIDS, NO_VAR)


