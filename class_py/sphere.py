import math

class Sphere:
	def __init__(self, r):
		self.radius = r
	def volume(self):
		return 4/3*math.pi*self.radius**3
	def surface_area(self):          #()裡放self
		return 4*math.pi*self.radius**2
s0 = Sphere(2)
print(s0.volume())
print(s0.surface_area())
