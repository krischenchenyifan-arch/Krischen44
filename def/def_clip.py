#函數內不只能放元素，串列也行
#把元素剪裁成介於minv與maxv間的數字
def clip(lst, minv, maxv):
	new_lst = []
	for n in lst:
		if n < minv:
			new_lst.append(minv)
		#else if n > maxv:（wrong）
		elif n > maxv:
			new_lst.append(maxv)
		else:
			new_lst.append(n)
	return new_lst
#主程式，測試函數
nums = [3, 7, 2, 9, 5, 12, 1]
print(clip(nums, 4, 10))
