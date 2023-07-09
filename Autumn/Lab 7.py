import numpy as np
import random
import math
import matplotlib.pyplot as plt

sigma_a = 0.03
H = 100
h = 1
D = 1
L = (1 / sigma_a) ** (1 / 2)
Q = 1
q = Q * h
a = - D / h
c = - D / h
N = int(H / h)

b = np.zeros(N)
b[0] = sigma_a * h - c
b[N - 1] = sigma_a * h - a + 0.5
b[1:N - 1] = sigma_a * h - c - a

A = np.zeros(N)
A[0] = q / b[0]

B = np.zeros(N)
B[0] = c / b[0]

for i in range(N - 1):
    A[i + 1] = (q - a * A[i]) / (b[i + 1] - a * B[i])
    B[i + 1] = c / (b[i + 1] - a * B[i])

Flow = np.zeros(N)
Flow[N - 1] = A[N - 1]
for i in range(N - 2, -1, -1):
    Flow[i] = A[i] - B[i] * Flow[i + 1]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
x = np.arange(H)
Flow_Theory = (Q / sigma_a) * (1 - np.cosh(x / L) / np.cosh(H / L))
ax.plot(x, Flow_Theory, color='blue', linewidth=1, label="Теория")
ax.plot(x, Flow, color='red', label="Модель")
plt.xlabel('Шаг пластины')
plt.ylabel('Поток')
plt.legend()
plt.grid()
plt.show()
