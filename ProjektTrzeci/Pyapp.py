import Matrix
import numpy as num

first = num.random.rand(50, 50)
second = num.random.rand(50, 50)
vec = num.random.rand(50)

res1, res2 = Matrix.gauss(first, vec, 50)

res3 = Matrix.multiplication(first, second, 50, 50)

print(res1)
print(res2)
print(res3)