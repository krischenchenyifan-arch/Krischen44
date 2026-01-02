import numpy as np
import random
def calc_fitness(x):
    return abs(x*x - 5)


population = np.array([0.2, 1.2, 2.0, 1.5])
# We should also store the fitness in an array
fitness = np.array([0.0, 0.0, 0.0, 0.0])


babies = np.array([0,0,0,0])
'''
population = {
    "A" : 0.2,
    "B" : 1.2,
    "C" : 2.0,
    "D" : 1.5
}
'''
ID = 0
for x in population:
    fitness[ID] = calc_fitness(x)
    ID += 1

# Now print out the fitness
print(fitness)

# Find the fittest using numpy's minimum function
fittest_ID  = np.argmin(fitness)
print(f"fittest_ID = {fittest_ID}")
print(f"Best solution: x = {population[fittest_ID]}, f(x) = {fitness[fittest_ID]}")

# Compute new baby X values
ID = 0
for x in population:
    Parent_X_ID = -1
    Parent_Y_ID = -1
    print(f"Computing baby value for {ID}")
    if (ID == fittest_ID):
        print("This is already the fittest; we need to randomly pick another parent")
        Parent_X_ID = ID
        # The other parent is randomly selected from the others
        #all_numbers = list(range(len(population)))
        #all_numbers.remove(Parent_X_ID)
        #Parent_Y_ID = random.choice(all_numbers)
        #result = 0.5 * (population[Parent_X_ID] + population[Parent_Y_ID])
        #print(f"{Parent_X_ID}, {Parent_Y_ID}")
        #print(f" x = {result}")
        while (Parent_Y_ID == -1):
            Test_ID = int(np.floor(np.random.rand()*4))
            if (Test_ID != ID):
                Parent_Y_ID = Test_ID

        print(f"{Parent_X_ID}, {Parent_Y_ID}")



    else:
        print("We need this parent, plus the best parent")
        Parent_X_ID = ID
        Parent_Y_ID = fittest_ID    
        print(f"{Parent_X_ID}, {Parent_Y_ID}")


    ID += 1