import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

class DifferentialEvolutionPeaks:
    def __init__(self, no_kids, no_var=2, no_gen=10, fr=1.0, sigma=1.0):
        self.NO_KIDS = no_kids
        self.NO_VAR = no_var
        self.NO_GEN = no_gen
        self.FR = fr
        self.SIGMA = sigma
        
        self.kid_dna = -3 + 6 * np.random.rand(self.NO_KIDS, self.NO_VAR)
        self.kid_fitness = np.zeros(self.NO_KIDS)
        self.history_fitness = np.empty(self.NO_GEN)
        self.best_index = -1

    def compute_fitness(self, x):
        X = x[0]
        Y = x[1]
        fitness = 3 * (1 - X)**2 * math.exp(-(X**2) - (Y + 1)**2) \
                - 10 * (X / 5 - X**3 - Y**5) * math.exp(-X**2 - Y**2) \
                - (1 / 3) * math.exp(-(X + 1)**2 - Y**2)
        return fitness

    def compute_best_kid(self, fitness_array):
        best_fitness = -10000.0
        best_index = -1
        for i in range(len(fitness_array)):
            if fitness_array[i] > best_fitness:
                best_fitness = fitness_array[i]
                best_index = i
        return best_index

    def max_mod(self, x):
        if x > 0.5:
            return x
        else:
            return 1.0 - x

    def compute_next_generation(self):
        trial_dna = np.empty(self.NO_VAR)
        new_dna = self.kid_dna.copy()
        new_fitness = self.kid_fitness.copy()

        for i in range(self.NO_KIDS):
            parent_a = self.best_index
            parent_b = parent_a
            parent_c = parent_a
            while (parent_a == parent_b) or (parent_a == parent_c) or (parent_b == parent_c):
                parent_b = np.random.randint(self.NO_KIDS)
                parent_c = np.random.randint(self.NO_KIDS)

            for j in range(self.NO_VAR):
                rf = self.FR * self.max_mod(np.random.rand())
                rnf = self.SIGMA * np.random.randn()
                trial_dna[j] = rf * self.kid_dna[parent_a, j] + (1.0 - rf) * self.kid_dna[parent_b, j] + rnf * (self.kid_dna[parent_b, j] - self.kid_dna[parent_c, j])

            trial_fitness = self.compute_fitness(trial_dna)

            if trial_fitness > self.kid_fitness[i]:
                new_fitness[i] = trial_fitness
                new_dna[i, :] = trial_dna

        self.kid_dna = new_dna
        self.kid_fitness = new_fitness

    def run(self):
        for i in range(self.NO_KIDS):
            self.kid_fitness[i] = self.compute_fitness(self.kid_dna[i, :])

        self.best_index = self.compute_best_kid(self.kid_fitness)

        for gen in range(self.NO_GEN):
            self.compute_next_generation()
            self.best_index = self.compute_best_kid(self.kid_fitness)
            self.history_fitness[gen] = self.kid_fitness[self.best_index]

        return self.history_fitness

# ==========================================
# 執行 1000 次模擬並計算平均適應值
# ==========================================
NO_GEN = 10
NUM_SIMULATIONS = 1000
population_sizes = [10, 100]

average_results_fitness = {}

fig, ax = plt.subplots(figsize=(10, 6))

for no_kids in population_sizes:
    print(f"正在執行 N = {no_kids} 的 {NUM_SIMULATIONS} 次模擬...")
    all_histories = np.zeros((NUM_SIMULATIONS, NO_GEN))
    
    for sim in range(NUM_SIMULATIONS):
        de = DifferentialEvolutionPeaks(no_kids=no_kids, no_gen=NO_GEN)
        all_histories[sim, :] = de.run()
        
    # 計算 1000 次模擬在每一個世代的平均值
    average_results_fitness[no_kids] = np.mean(all_histories, axis=0)

# 繪製平均適應值對世代圖
for no_kids in population_sizes:
    ax.plot(average_results_fitness[no_kids], marker='o', label=f'Average N = {no_kids}')

ax.set(xlabel='Generation Number', ylabel='Average Best Fitness', title=f'Average Fitness vs Generation ({NUM_SIMULATIONS} Simulations)')
ax.legend()
ax.grid(True)
plt.show()
