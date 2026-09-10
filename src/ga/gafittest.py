import random

# --- 函數設定 ---

# (A) Fitness 函數：計算目前的 x 與 目標值(sqrt(5)) 的接近程度
def fitness(x):
    # 流程圖中為 abs(x * x - 5)
    return abs(x * x - 5)

# (B) CalcFittest 函數：從群體中找出適應度最好（數值最小）的索引
def CalcFittest(F):
    bestF = 10000
    bestID = -1
    # NO_P 是群體數量
    for p in range(len(F)):
        if F[p] < bestF:
            bestF = F[p]
            bestID = p
    return bestID

# --- 主程式開始 ---

# 初始化參數
G = 0            # 目前世代
NO_G = 10        # 最大世代數 (世代數限制)
NO_P = 4         # 群體大小 (Population size)

# 初始化群體 x 和其對應的適應度 F
# 流程圖標示為 [rand, rand, rand, rand]
x = [random.uniform(0, 5) for _ in range(NO_P)]
F = [fitness(val) for val in x]

# 找出初始群體中最優者
Fittest = CalcFittest(F)


while G < NO_G:
#根據條件判斷，決定是否重複或停止，用法為「while 條件:」
#如果條件判斷為 True，就會不斷執行迴圈內容，如果判斷為 False，就會停止迴圈
    new_x = [0] * NO_P
    new_F = [0] * NO_P
    
    
    for p in range(NO_P):
        #for 變數 in 可迭代的物件，for 迴圈會依序將可以迭代的物件取出，賦值給指定的變數
        # 選擇父母 (Parent Selection)
        if p == Fittest:
            
            parentA = p
            choices = [i for i in range(NO_P) if i != p]
            parentB = random.choice(choices)
        else:
            
            parentA = p
            parentB = Fittest
            
        
        # 公式：x* = 0.5 * (x[A] + x[B]) + 0.5 * random
        x_star = 0.5 * (x[parentA] + x[parentB]) + 0.5 * (random.random() * 2 - 1)
        f_star = fitness(x_star)
        
        
        if f_star < F[p]:
            new_x[p] = x_star
            new_F[p] = f_star
        else:
            new_x[p] = x[p]
            new_F[p] = F[p]
            
    
    x = new_x
    F = new_F
    G += 1
    Fittest = CalcFittest(F)
    
    
    #print(f"Generation {G}: Best X = {x[Fittest]:.4f}, Fitness = {F[Fittest]:.4f}")

#print("-" * 30)
#print("Final Results:")
print("x values:", x)
print("F values:", F)
#print(f"Best Solution x: {x[Fittest]}")