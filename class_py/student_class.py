class Student:
	student_count = 0	#類別屬性
	def __init__(self, name):
		self.name = name
		Student.student_count += 1
		
	@classmethod
	def get_student_count(cls):
		return f"There are {cls.student_count} students."
s1 = Student('Alice')
s2 = Student('bob')
s3 = Student('Charlie')
print(Student.get_student_count())
