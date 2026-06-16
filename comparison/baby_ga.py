# Baby_ga.py
# This is the introduction to Genetic Algorithms
# with very simple reproduction and no mutation.

import numpy as np
import tools_ga as tools

def ComputeNextGen_GA1(DNA, FITNESS, BESTINDEX):
    NO_KIDS, NO_VAR = DNA.shape
    trial_dna = np.empty(NO_VAR)
    new_dna = DNA.copy()
    new_fitness = FITNESS.copy()
 
    for i in range(NO_KIDS):
        Parent_A = BESTINDEX
        Parent_B = Parent_A
        while (Parent_A == Parent_B):
            Parent_B = np.random.randint(NO_KIDS)

        # X_baby = 1/2 * (X_A + X_B)
        for j in range(NO_VAR):
            W = 0.5
            trial_dna[j] = W * DNA[Parent_A, j] + (1.0 - W) * DNA[Parent_B, j]

        trial_fitness = tools.ComputeFitness(trial_dna)
        if (np.absolute(trial_fitness) < np.absolute(FITNESS[i])):
            new_fitness[i] = trial_fitness
            new_dna[i, :] = trial_dna
    return new_dna, new_fitness


