import math

def isprime(x):
    # 質數必須大於 1
    if x <= 1:
        return False
    
    # 2 是唯一的偶數質數
    if x == 2:
        return True
    
    # 排除 2 以外的所有偶數
    if x % 2 == 0:
        return False
    
    # 從 3 開始，只檢查奇數因數，一路檢查到根號 x（包含根號 x 本身）
    # 使用 math.isqrt(x) 可以直接得到整數平方根，效率更高
    for i in range(3, math.isqrt(x) + 1, 2):
        if x % i == 0:
            return False  # 如果能被整除，就不是質數
            
    return True  # 順利通過檢查，代表是質數

# ==========================================
# 實際測試
# ==========================================
test_numbers = [0, 2, 5, 12, 23, 37, 107]

for num in test_numbers:
    print(f"數字 {num:2d} 是否為質數？ {isprime(num)}")
