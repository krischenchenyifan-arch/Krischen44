class Student:
	count = 0 #屬於類別屬性，用於計數，由所有物共享
	def __init__(self, name , age):
		self.name = name
		self.age = age
		Student.count += 1   #每次建立物件時計數加一
	@classmethod
	def from_age(cls, name, age):
		return cls(name, age)
	@classmethod
	def get_count(cls):
		return cls.count

print(f"total :{Student.get_count()}")
s1 = Student('Alice', 20)
s2 = Student.from_age('Bob', 22)
print(f"Name: {s1.name}, Age: {s1.age}")
print(f"Name: {s2.name}, Age: {s2.age}")
print(f"total: {Student.get_count()}")
print(f"total: {s1.get_count()}")
