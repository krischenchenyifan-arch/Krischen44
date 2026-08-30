class Student:
	def __init__(self, name, email):
		self.name = name
		self.email = email
	@staticmethod
	def check(email):
		return '正確的email' if '@' in email else '錯誤的email'

print(Student.check('invalid-email'))
#class呼叫
s1 = Student('Alice', 'alice@email.com')
print(Student.check(s1.email))
#s1.email是self.email
#instance呼叫
print(s1.check('alice@email.com'))
