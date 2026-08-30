class Temperature:
    def __init__(self, fahrenheit):
        self.fahrenheit = fahrenheit  # 儲存華氏溫度到實體屬性

    @classmethod
    def from_celsius(cls, celsius):
        # 1. 先把攝氏換算成華氏
        fahrenheit = (9 / 5) * celsius + 32
        # 2. 關鍵：用 cls(值) 來建立並回傳一個 Temperature 實體物件
        return cls(fahrenheit)

# 主程式
temp = Temperature.from_celsius(25)
print(temp.fahrenheit)
