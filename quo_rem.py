print(f"======Function to calculate quotient(商) and remainder(餘數)======")

def quo_rem(a, b):
        quo = a // b
        #計算商
        rem = a % b
        #計算餘數
        return quo, rem
        #傳回商和餘數

a, b = 100000, 93
q, r = quo_rem(a, b)
print(f"{a}={b}*{q}+{r}")
