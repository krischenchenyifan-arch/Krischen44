class PrimeChecker:
	@staticmethod
	def is_prime(n):
		for i in range(2, int(n**0.5) + 1):
			if n % i == 0:
				return False
		return True
	@classmethod
	def find_primes(cls, limit):
		primes = [n for n in range(2, limit + 1) if cls.is_prime(n) ]
		return primes
	@classmethod
	def prime_append(cls, limit):
		primes = []
		for n in range(2, limit + 1):
			if cls.is_prime(n):
				primes.append(n)
		return primes

print(PrimeChecker.is_prime(23))
print(PrimeChecker.is_prime(47))
limit = 100
primes = PrimeChecker.find_primes(limit)
#串列式生成(List Comprehension)
print(f'小於等於{limit}的質數有：{primes}')
print(f"小於等於{limit}的質數有：{PrimeChecker.prime_append(limit)}")
