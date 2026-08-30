import random

class Pet:
	def __init__(self, name):
		self.name = name
	def race(self, other):
		my_speed = random.randint(1, 10)
		other_speed = random.randint(1, 10)
		winner = self.name if my_speed > other_speed else other.name
		#print(f"{winner}跑贏了～")
		print(f"小橘貓:{my_speed}")
		print(f"虎斑貓:{other_speed}")
		return f"{winner}跑贏了～"
pet1 = Pet("小橘貓")
pet2 = Pet("虎斑貓")
result = pet1.race(pet2)
print(result)
