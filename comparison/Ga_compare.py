import matplotlib.pyplot as plt
import numpy as np
import math

import tools_ga as tools
import baby_ga as baby

# ---- 演算法 2：GA (#2) 三父母含變異機制 ----
def ComputeNextGen_GA2(DNA, FITNESS, BESTINDEX, FR, SIGMA):
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


# ==================== 主程式開始 ====================

# 基礎參數設定
NO_KIDS = 30
NO_VAR = 1   # 一維問題
NO_GEN = 15  # 稍微拉長代數到 15 代，觀察收斂對比
FR = 1.0
SIGMA = 0.5  # 調整擾動標準差，避免一維問題跳躍過大

# 1. 產生初始群體 DNA (起跑點)
initial_dna = 5.0 * np.random.rand(NO_KIDS, NO_VAR)
initial_fitness = np.zeros(NO_KIDS)
for i in range(NO_KIDS):
    initial_fitness[i] = tools.ComputeFitness(initial_dna[i, :])

# ---- 分流複製給 GA1 與 GA2 ----
dna_ga1 = initial_dna.copy()
fit_ga1 = initial_fitness.copy()
best_ga1 = tools.ComputeBestKid(fit_ga1)

dna_ga2 = initial_dna.copy()
fit_ga2 = initial_fitness.copy()
best_ga2 = tools.ComputeBestKid(fit_ga2)

# 歷史紀錄陣列
history_fit_ga1 = np.empty(NO_GEN)
history_dna_ga1 = np.zeros((NO_GEN, NO_VAR))

history_fit_ga2 = np.empty(NO_GEN)
history_dna_ga2 = np.zeros((NO_GEN, NO_VAR))

# 2. 開始同時演化
for gen in range(NO_GEN):
    # 執行 GA (#1)
    dna_ga1, fit_ga1 = baby.ComputeNextGen_GA1(dna_ga1, fit_ga1, best_ga1)
    best_ga1 = tools.ComputeBestKid(fit_ga1)
    history_fit_ga1[gen] = fit_ga1[best_ga1]
    history_dna_ga1[gen, :] = dna_ga1[best_ga1, :]
    
    # 執行 GA (#2)
    dna_ga2, fit_ga2 = ComputeNextGen_GA2(dna_ga2, fit_ga2, best_ga2, FR, SIGMA)
    best_ga2 = tools.ComputeBestKid(fit_ga2)
    history_fit_ga2[gen] = fit_ga2[best_ga2]
    history_dna_ga2[gen, :] = dna_ga2[best_ga2, :]


# ==================== 繪圖呈現 ====================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# 圖左：DNA 估計值的收斂路徑比較
ax1.plot(history_dna_ga1, 'b-o', label='GA (#1) - 2 Parents')
ax1.plot(history_dna_ga2, 'r-s', label='GA (#2) - 3 Parents + Var')
ax1.axhline(y=math.sqrt(5), color='g', linestyle='--', label='Target sqrt(5)')
ax1.set_xlabel('Generation Number')
ax1.set_ylabel('Estimated DNA Value (x)')
ax1.set_title('DNA Value (x) Convergence')
ax1.legend()
ax1.grid(True)

# 圖右：適應度絕對值差距下降曲線比較 (|Fitness|)
ax2.plot(np.absolute(history_fit_ga1), 'b-o', label='GA (#1) - 2 Parents')
ax2.plot(np.absolute(history_fit_ga2), 'r-s', label='GA (#2) - 3 Parents + Var')
ax2.set_yscale('log')  # 使用對數座標，更能看出微小差距的逼近
ax2.set_xlabel('Generation Number')
ax2.set_ylabel('Absolute Fitness |5 - x^2|')
ax2.set_title('Fitness Error Decay (Log Scale)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('comparision.png')


# ==================== INDIVIDUAL FINAL REPORTS ====================

print("\n" + "="*20 + " FINAL REPORT: GA (#1) " + "="*20)
print("Average fitness of final generation = %g" % np.mean(fit_ga1))
print("Best Kid's Fitness (5 - x^2)        = %g" % fit_ga1[best_ga1])
print("Best Kid's DNA (Estimated sqrt(5))  = %g" % dna_ga1[best_ga1, 0])
print("Error compared to real sqrt(5)      = %g" % np.abs(dna_ga1[best_ga1, 0] - math.sqrt(5)))

print("\n" + "="*20 + " FINAL REPORT: GA (#2) " + "="*20)
print("Average fitness of final generation = %g" % np.mean(fit_ga2))
print("Best Kid's Fitness (5 - x^2)        = %g" % fit_ga2[best_ga2])
print("Best Kid's DNA (Estimated sqrt(5))  = %g" % dna_ga2[best_ga2, 0])
print("Error compared to real sqrt(5)      = %g" % np.abs(dna_ga2[best_ga2, 0] - math.sqrt(5)))

print("\n" + "="*23 + " BENCHMARK TARGET " + "="*23)
print("Real sqrt(5) value                 = %g" % math.sqrt(5))
print("="*64)
