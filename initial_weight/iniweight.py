import numpy as np
#類別的基本程式練習 ——create random weight
class Weight:
	#初始化特徵(Attributes)
	#本情況中有self和id兩種特徵
	def __init__(self, id):
		self.value = np.random.rand()
		self.id = id
	#Method:物件裡的行為，也就是class裡面的函式
	def show_details(self):
		print(f"id = {self.id}, value = {self.value}")

weights = []
NO_weight = 16
for i in range(NO_weight):
	weight = Weight(i)
	weights.append(weight)

for weight in weights:
	Weight.show_details(weight)

class Cat:
	def __init__(self, name, color):
		self.name = name
		self.color = color
	def meow(self):
		print(f"{self.name} say: meow~ I am {self.color} cat!")
		#return f"{self.name} say: meow~ I am {self.color} cat!"
cat1 = Cat("Dave", "yellow")
cat2 = Cat("Chen", "red")
print(cat1.name)
print(cat2.color)
#print(cat1.meow())
cat1.meow()
#如果函式裡有print了，外面不用print
#函式裡用return，外面要print
