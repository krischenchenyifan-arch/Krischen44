import copy

'''
#1.dict practice
student = {
	'name' : 'Bob',
	'age' : 31,
	'job' : 'programer'
}
print(f"Keys for student are {student.keys()}")
print(f"The type is {type(student)}")
print(f"The name is {student['name']}")
student['salary'] = 3000
print(f"Object is : {student}")

'''


'''
#2.list基本操作
grade = [[67,80,87,69],
	 [71,80,65,53],
	 [77,58,60,49]]
print(len(grade))   #查詢grade的長度，得到3，代表grade 中有3個元素
print(len(grade[0]))
print(grade[0])
a = [c[0] for c in grade]
print(a)
'''

'''
#計算學生平均成績
students =[
	['Alice', [85,90,78]],
	['Bob', [88,76,92]],
	['Charlie', [90,85,89]]
]
for student in students:
	name = student[0]
	No_score = len(student[1])
#	average_score = sum(c for c in student[1])/ No_score
	average_score = sum(student[1]) / No_score
	print(f"{name} : {average_score:.2f}")
#使用sum()計算成績總和
#也可以寫{(average_score):.2f}
#:.2f代表取數值至小數點後兩位

'''

'''
#copy
lst = [[1,2], 34]
#1.設定運算：(指向同一個lst)
a = lst
a[0][0] = 99
print('設定運算：', lst)

#2.淺拷貝.copy()

lst = [[1,2], 34]
b = lst.copy()
b[0][0] = 99
print('淺拷貝：', lst)
#輸出:[[99,2], 34]，lst仍有受到影響

#3.深拷貝copy.deepcopy()

lst = [[1,2], 34]
c = copy.deepcopy(lst)
#須先載入copy模組：import copy
c[0][0] = 99
print('深拷貝：', lst)
#輸出:[[1,2], 34]，lst不受到影響

'''

a = round(5.25 , 1)
b = round(5.35 , 1)   #5.3??????
c = round(3.35 , 1)   #3.4
d = pow(2, 6)
e = min(4, 5, 3)
f = abs(-4.8)
g = max(4, 5, 3)
h = round(1.35, 1)    #1.4
print(a, b, c, d, e, f, g, h)
 
