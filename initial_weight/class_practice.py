#calculate circle area, perimeter and volume using fuction and class separately.
#calculate with function.
import math
radius =4
def area(radius):
	return math.pi*radius**2 
def perimeter(radius):
	return 2*math.pi*radius
def volume(radius):
	return 4/3*math.pi*radius**3

for radius in range(4,11):
	c_area = area(radius)
	c_perimeter = perimeter(radius)
	c_volume = volume(radius)

	print(f"半徑：{radius},面積： {c_area},週長： {c_perimeter},體積：{c_volume} ")
print(f"=====using class=====")
class Circle:
	def __init__(self, radius2=4):
		self.radius = radius2
	def area2(self):
		return math.pi*(self.radius)**2
	def perimeter2(self):
		return 2*math.pi*self.radius
	def volume2(self):
		return 4/3*math.pi*(self.radius)**3

for radius2 in range(4,11):
	c1 = Circle(radius2)
	print(f"半徑：{c1.radius},面積： {c1.area2()},週長： {c1.perimeter2()},體積：{c1.volume2()}")

#add ()after the function in order to run the method
#c1.area2(), c1.perimeter2(), c1.volume2()
