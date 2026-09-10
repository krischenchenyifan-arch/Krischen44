#傳回1+2+3+4+...+n
def add(n):
	num = 0
	for i in range(n+1):
		num += i
	return num
print(add(200))
