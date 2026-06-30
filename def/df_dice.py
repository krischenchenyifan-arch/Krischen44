import random
def roll_dice():
#函數名稱
#沒有參數
	num = random.randint(1, 6)
	print(f"The number is {num}")
	#return f"The number is {num}"
	#函數本體
	#沒有回傳值no return
roll_dice()
#print(roll_dice())
#主程式
#呼叫函數
print(f"=======增加次數=======")
def roll_dice_n(n):
	for _ in range(n):
		num2 = random.randint(1, 6)
		print(f"The number is {num2}")

roll_dice_n(3)

print(f"======傳出總點數======")

def roll_dice_total(n):
	total = 0
	for _ in range(n):
		result = random.randint(1, 6)
		print(f"The number is {result}")
		total += result
		#total = total + result
	return total
	#傳回總點數
n = 3
total_dice = roll_dice_total(n)
print(f"{n}次總和為 {total_dice}")

