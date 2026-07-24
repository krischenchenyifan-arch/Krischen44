import math
class NumberUtils:
	def is_prime(n):
		for i in range(2, math.floor(n**0.5) + 1):
			if n % i == 0:
				return False
		return True
print(NumberUtils.is_prime(37))
print(NumberUtils.is_prime(45))


'''
class Primechecker2:
	@staticmethod
	def is_prime(n):
		for i in range(2, int(n**0.5) + 1):
			if n % i == 0:
				return False
		return True
	@classmethod
	def find_prime(cls, limit):
		primes = [n for n in range(2, limit +1) if cls.is_prime(n)]
		return primes

limit = 20
print(Primechecker2.find_prime(limit))
How to get a integer?
int(x)		無條件捨去
round(x)	四捨五入
math.floor(x)	無條件捨去
math.ceil(x)	無條件進位
'''
