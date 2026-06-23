import numpy as np
class Weight:
	def __init__(self, id):
		self.value = np.random.rand()
		self.id = id
	def show_details(self):
		print(f"id = {self.id}, value = {self.value}")

weights = []
NO_weight = 16
for i in range(NO_weight):
	weight = Weight(i)
	weights.append(weight)

for weight in weights:
	Weight.show_details(weight)
