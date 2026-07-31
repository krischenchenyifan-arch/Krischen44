import math

input = [5.1, 3.8, 1.6, 0.2]

bias = {
	0: 0, 1: 0, 2: 0, 3: 0,
	4: -0.0260747, 5: 0.339481, 6: 0.542879, 7: -0.729046, 8: 0.806341,
	9: 0.4566, 10: 0.475939, 11: -0.564546
}

weight = {
	0: 7.52436, 1: 8.84869, 2: -19.8074, 3: -3.54702, 4: 5.22944,
	5: -2.48485, 6: 6.53303, 7: -14.6117, 8: -22.8034, 9: -4.34393,
	10: 5.99496, 11: 0.207013, 12: -1.11588, 13: -6.96244, 14: -18.3058,
	15: -5.75681, 16: -5.80779, 17: -0.535812, 18: -7.00339, 19: -8.97093,
	20: 12.2787, 21: -7.31267, 22: 1.44651,
	23: -4.53042, 24: 7.34816, 25: -5.36046,
	26: -0.124184, 27: 5.96869, 28: -5.73783,
	29: -2.75999, 30: -12.3828, 31: 0.114309,
	32: 4.90098, 33: -10.149, 34: -8.4269
}
def sigmoid(z):
	return 1/(1+math.exp(-z))

class Neuron:
	def __init__(self, input, weight):
		self.input = input
		self.weight = weight
	def output(self):
		z = self.input*self.weight
		return z


z4_0 = Neuron(input[0], weight[0])
z4_1 = Neuron(input[1], weight[5])
z4_2 = Neuron(input[2], weight[10])
z4_3 = Neuron(input[3], weight[15])
z4 = z4_0.output() + z4_1.output() + z4_2.output() + z4_3.output() + bias[4]
z4_final = sigmoid(z4)
print(f"Neuron 4 output is {z4_final}")

z5_0 = Neuron(input[0], weight[1])
z5_1 = Neuron(input[1], weight[6])
z5_2 = Neuron(input[2], weight[11])
z5_3 = Neuron(input[3], weight[16])
z5 = z5_0.output() + z5_1.output() + z5_2.output() + z5_3.output() + bias[5]
z5_final = sigmoid(z5)
print(f"Neuron 5 output is {z5_final}")

z6_0 = Neuron(input[0], weight[2])
z6_1 = Neuron(input[1], weight[7])
z6_2 = Neuron(input[2], weight[12])
z6_3 = Neuron(input[3], weight[17])
z6 = z6_0.output() + z6_1.output() + z6_2.output() + z6_3.output() + bias[6]
z6_final = sigmoid(z6)
print(f"Neuron 6 output is {z6_final}")

z7_0 = Neuron(input[0], weight[3])
z7_1 = Neuron(input[1], weight[8])
z7_2 = Neuron(input[2], weight[13])
z7_3 = Neuron(input[3], weight[18])
z7 = z7_0.output() + z7_1.output() + z7_2.output() + z7_3.output() + bias[7]
z7_final = sigmoid(z7)
print(f"Neuron 7 output is {z7_final}")

z8_0 = Neuron(input[0], weight[4])
z8_1 = Neuron(input[1], weight[9])
z8_2 = Neuron(input[2], weight[14])
z8_3 = Neuron(input[3], weight[19])
z8 = z8_0.output() + z8_1.output() + z8_2.output() + z8_3.output() + bias[8]
z8_final = sigmoid(z8)
print(f"Neuron 8 output is {z8_final}")

print(f"======start calculate output layer ouputs======")

z9_4 = Neuron(z4_final, weight[20])
z9_5 = Neuron(z5_final, weight[23])
z9_6 = Neuron(z6_final, weight[26])
z9_7 = Neuron(z7_final, weight[29])
z9_8 = Neuron(z8_final, weight[32])
z9 = z9_4.output() + z9_5.output() + z9_6.output() + z9_7.output() + z9_8.output() + bias[9]
print(f"Neuron 9 output is {sigmoid(z9)}")

z10_4 = Neuron(z4_final, weight[21])
z10_5 = Neuron(z5_final, weight[24])
z10_6 = Neuron(z6_final, weight[27])
z10_7 = Neuron(z7_final, weight[30])
z10_8 = Neuron(z8_final, weight[33])
z10 = z10_4.output() + z10_5.output() + z10_6.output() + z10_7.output() + z10_8.output() + bias[10]
print(f"Neuron 10 output is {sigmoid(z10)}")

z11_4 = Neuron(z4_final, weight[22])
z11_5 = Neuron(z5_final, weight[25])
z11_6 = Neuron(z6_final, weight[28])
z11_7 = Neuron(z7_final, weight[31])
z11_8 = Neuron(z8_final, weight[34])
z11 = z11_4.output() + z11_5.output() + z11_6.output() + z11_7.output() + z11_8.output() + bias[11]
print(f"Neuron 11 output is {sigmoid(z11)}")
