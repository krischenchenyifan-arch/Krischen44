class Circle:
	pi = 3.14159
	#pi是類別屬性(class atribute)，放在def __init__():外面
	def __init__(self, radius):	#radius是實例屬性(instance attribute)
		self.radius = radius
	def area(self):
		return Circle.pi*self.radius**2
	def perimeter(self):
		return 2*Circle.pi*self.radius
c1 = Circle(5)
print(f"area:{c1.area():.2f}")
print(f"perimeter:{c1.perimeter():.2f}")
Circle.pi = 3.14
print(f"area:{c1.area():.2f}")
print(f"perimeter:{c1.perimeter():.2f}")
