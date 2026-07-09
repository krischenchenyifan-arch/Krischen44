class Employee:
	company = 'Unknown'
	def __init__(self, name, salary):
		self.name = name
		self.salary = salary
	def display_info(self):
	#def display_info(self, company):
		print(f"name:{self.name}, salary:{self.salary}, company:{Employee.company}")

Employee.company = 'TechCorp'
emp1 = Employee('Alice', 50000)
emp2 = Employee('Bob', 60000)
emp1.display_info()
emp2.display_info()
#emp1.display_info('TechCorp')
#emp2.display_info('TechCorp')
#TypeError: Employee.display_info() missing 1 required positional argument: 'company'
#為何刪掉Employee.company = 'TechCorp'後輸出仍然是company = 'Unknown'
#??????
