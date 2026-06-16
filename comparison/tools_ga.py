# tools_ga.py
# Contains things used by all GA solvers,
# like fitness calculation.

def ComputeFitness(x):
    # 目標函數：f(x) = 5 - x^2
    X = x[0]
    return 5.0 - X**2
