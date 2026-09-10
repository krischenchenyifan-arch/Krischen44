class Circle:
	def __init__(self, radius = 5, pi = 3.14):
		self.radius = radius
		self.pi = pi
	def area(self):
		return self.pi*self.radius**2
	def perimeter(self):
		return 2*self.pi*self.radius
#print(f"{Circle.area():.2f}")
#print(f"{Circle.perimeter():.2f}")
a1 = Circle(6, 3.14159)
print(f"{a1.area():.2f}")
print(f"{a1.perimeter():.2f}")
