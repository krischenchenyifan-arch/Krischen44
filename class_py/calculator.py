import math

class Calculator:
	def __init__(self, a, b):
		self.n1 = a
		self.n2 = b
	def add(self):
		return self.n1 + self.n2
	def gcd(self):
		return math.gcd(self.n1, self.n2)
	def lcm(self):
		return math.lcm(self.n1, self.n2)
	def power(self):
		#pow(base, exp)
		return pow(self.n1, self.n2)

c0 = Calculator(2, 10)
print(c0.add())
print(c0.gcd())
print(c0.lcm())
print(c0.power())
