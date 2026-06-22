# tools_ga.py
# Contains things used by all GA solvers,
# like fitness calculation.

import numpy as np

def MaxMod(X):
    # GA (#2) 需要使用的權重調節函數
    if (X > 0.5):
        return X
    else:
        return (1.0 - X)

def ComputeFitness(x):
    # 目標函數：f(x) = 5 - x^2
    X = x[0]
    return 5.0 - X**2

def ComputeBestKid(FITNESS):
    # 尋找適應度絕對值最接近 0 (最優秀) 的個體索引
    NO_KIDS = len(FITNESS)
    BestFitness = 10000.0
    BestIndex = -1
    for i in range(NO_KIDS):
        if (np.absolute(FITNESS[i]) < BestFitness):
            BestFitness = np.absolute(FITNESS[i])
            BestIndex = i
    return BestIndex
