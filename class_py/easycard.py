class Easycard:
	def __init__(self, owner, amount = 0 ):
		self.name = owner
		self.balance = amount
		print(f"{self.name}儲值{amount}")
	def add_value(self, amount):
		self.balance += amount
		print(f"{self.name}儲值{amount}，餘額：{self.balance}")
	def spend(self, amount):
		if amount <= self.balance:
			self.balance -= amount
			print(f"{self.name}消費{amount}，餘額：{self.balance}")
		else:
			print("餘額不足")
ac1 = Easycard("Mary", 500)
ac1.add_value(5000)
ac1.spend(80000)

ac2 = Easycard("Chen", 0)
ac2.spend(400)
ac2.add_value(600)
ac2.spend(400)


