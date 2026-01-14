import random

def fitness(x):
    return abs(x*x - 5)

def calcFittest(F):
    bestF = 10000
    bestID = -1
    for p in range(NO_P):
        if F[p] < bestF:
            bestF = F[p]
            bestID = p
    return bestID

G = 0
NO_P = 4
NO_G = 10


#x[NO_P] = [rand, rand, rand, rand]
#F[NO_P] = fitness(x)
fittest = calcFittest(F)
#print(f"Generation {G}: Best X = {x[Fittest]:.4f}, Fitness = {F[Fittest]:.4f}")

#print("-" * 30)
#print("Final Results:")
#print(f"Best Solution x: {x[Fittest]}")