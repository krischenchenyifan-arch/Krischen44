import numpy as np



arr = np.array([
    [9, 5, 6, 8],
    [2, 4, 7, 6],
    [4, 5, 7, 4]
])
#矩陣需要用[]括起來方能印出

print(arr)                      #[[9, 5, 6, 8]
                                # [2, 4, 7, 6]
                                # [4, 5, 7, 4]]
print("=====Hello, World!=====")#=====Hello, World!=====
name = "Alice"                  #Alice
print(name)                     #print variables
print(123)                      #123     numbers
print(1 + 3)                    #4       expressions
print(2 * 3)                    #6       expressions
age = 30
print("Name :", name,           #輸出：  Name : Alice Age : 30     如何換行？
      "Age :", age)
print("Name =", name, "\nAge :", age)    #用\n (常用)
print("Name =", name)                    #分兩個print
print("Age =", age)
print("Mame =", name, sep="")            #用sep參數
print("Age =", age, sep="")
print("Name =", name, "\nAge =", age, sep="")
print(f"""Name : {name}                          
Age : {age}""")                          #f-string
#np.array = np.arange(15).reshape(3, 5)
#np.array([[ 9, 5, 8, 9],
          #[ 2, 3, 4, 7],
          #[ 3, 6, 9, 1]])

#print(f"{np.array}")
