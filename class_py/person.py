class Person:
	def __init__(self, name, age):
		self.name = name
		self.age = age
	def greet(self):
		#return f"Hi,I am {self.name},{self.age}"
		print(f"Hi,I am {self.name},{self.age}")
p = Person('Alice', 21)
#print(p.greet())
p.greet()
#有return 就要print
