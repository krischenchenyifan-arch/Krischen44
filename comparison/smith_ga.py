import numpy as np
import tools_ga as tools

def ComputeNextGen_GA(DNA, FITNESS, BESTINDEX, FR, SIGMA):
    NO_KIDS, NO_VAR = DNA.shape
    trial_dna = np.empty(NO_VAR)
    new_dna = DNA.copy()
    new_fitness = FITNESS.copy()
    
    for i in range(NO_KIDS):
        Parent_A = BESTINDEX
        Parent_B = Parent_A
        Parent_C = Parent_A
        while ((Parent_A == Parent_B) or (Parent_A == Parent_C) or (Parent_B == Parent_C)):
            Parent_B = np.random.randint(NO_KIDS)
            Parent_C = np.random.randint(NO_KIDS)

        # X_baby = X_baby + randn() * (X_B - X_C) 邏輯
        for j in range(NO_VAR):
            Rf = FR * tools.MaxMod(np.random.rand())
            Rnf = SIGMA * np.random.randn()
            trial_dna[j] = Rf * DNA[Parent_A, j] + (1.0 - Rf) * DNA[Parent_B, j] + Rnf * (DNA[Parent_B, j] - DNA[Parent_C, j])

        trial_fitness = tools.ComputeFitness(trial_dna)
        if (np.absolute(trial_fitness) < np.absolute(FITNESS[i])):
            new_fitness[i] = trial_fitness
            new_dna[i, :] = trial_dna
    return new_dna, new_fitness
