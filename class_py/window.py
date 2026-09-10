class Window:
	def __init__(self, width = 10, height = 5):
		self.width = width
		self.height = height
	def area(self):
		area = self.width*self.height
		print(area)
		#return area
w0 = Window()
w1 = Window(12, 8)
w0.area()
w1.area()
#print(w0.area())
#print(w1.area())

