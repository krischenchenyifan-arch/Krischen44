import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from GA import ComputeNextGeneration_Storn, ComputeNextGeneration_Smith, ComputeBestKid


#Define 2D Ackley's function(2 variances)
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
#[(-5<X<5), (-5<Y<5)] -------> kid_dna
#Global minimum exists at (0, 0) with f(X,Y) = 0
#====================================
#hyperparameters(超參數)
#(MF, CR)---> Storn , (FR, SIGMA)---> Smith
#====================================
#畫出DE演算法收斂曲線
#Convergence towards the global minima for Ackley's function 
solvers = [
	"Smith", "Smith", "Smith", "Smith",
	"Storn", "Storn", "Storn", "Storn"
]

parameter1 = [1.0, 1.0, 1.5, 1.5, 0.25, 0.25, 0.5, 0.5] #Smith:FR   \ Storn:MF
parameter2 = [1.0, 1.5, 1.0, 1.5, 0.8, 1.0, 0.8, 1.0]   #      SIGMA\	    CR

labels = [
	"Proposed DE, F = 1.0, \u03c3 = 1.0",
	"Proposed DE, F = 1.0, \u03c3 = 1.5",
	"Proposed DE, F = 1.5, \u03c3 = 1.0",
	"Proposed DE, F = 1.5, \u03c3 = 1.5",
	"DE, F = 0.25, CR = 0.8",
	"DE, F = 0.25, CR = 1.0",
	"DE, F = 0.5, CR = 0.8",
	"DE, F = 0.5, CR = 1.0"
]
fig, ax = plt.subplots(figsize = (10, 8))

for i in range(len(solvers)):
	current_solver = solvers[i] 
	p1 = parameter1[i]
	p2 = parameter2[i]
	current_label = labels[i]
	for run in range(NUM_RUNS):
		kid_dna = -5 + 10 * np.random.rand(NO_KIDS, NO_VAR)
	        #-5開始，長度10
		kid_fitness = np.empty(NO_KIDS)
		for j in range(NO_KIDS):
			kid_fitness[j] = Computefitness_Ackley(kid_dna[j,:])
		bestkid = ComputeBestKid(kid_fitness)
		
		for gen in range(NO_GEN):
			if (current_solver == "Smith"):
				kid_dna, kid_fitness = ComputeNextGeneration_Smith(kid_dna, kid_fitness, p1, p2, Computefitness_Ackley, None)
			elif (current_solver == "Storn"):
				kid_dna, kid_fitness = ComputeNextGeneration_Storn(kid_dna, kid_fitness, p1, p2, Computefitness_Ackley, None)
			

