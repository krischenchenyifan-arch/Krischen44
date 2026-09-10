class Factor:
	factor_list = [2, 3, 6, 8]
	factor_list1 = []
	def __init__(self, num):
		print(f'初始的factor_list: ',Factor.factor_list)
		self.num = num
	def find_factors(self):
		for i in Factor.factor_list:
			if self.num % i ==0:
				Factor.factor_list1.append(i)
		return Factor.factor_list1
	@classmethod
	def add_factors(cls, lst = []):
		for i in cls.lst:
			if 
		return cls()
	@classmethod
	def remove_factors(cls, lst):
		

	@classmethod
	def show_factor_list(cls):
	

	@staticmethod
	def isfactor(num, n):

f0 = Factor(12)
print(f0.find_factors())
Factor.add_factors([4, 6, 9])
Factor.show_factor_list()
Factor.remove_factors([3, 8])
Factor.show_factor_list()
print(Factor.isfactor(15, 3))
print(Factor.isfactor(15, 4))
