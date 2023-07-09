import numpy as np
import random
import math
import matplotlib.pyplot as plt

sigma_a = 0.03
sigma_gen = 0.0303
nu_f = 2.4
sigma_f = sigma_gen / nu_f
sigma_s = 0.3
sigma_c = sigma_a - sigma_f
H = 100
h = 1
D = 1
L = (1 / sigma_a) ** (1 / 2)
a = - D / h
c = - D / h
N = int(H / h)
Q = np.zeros(N)
Flow = np.zeros(N)
eps = 10e-4
Flow[0:N] = 1
q = np.zeros(N)
k_eff = 100
k = 1
while not abs(k - k_eff) < eps:
    Last_Norm = np.linalg.norm(Flow)
    k = k_eff
    Q[0:N] = sigma_f * nu_f * Flow[0:N] / k_eff
    q[0:N] = Q[0:N] * h
    b = np.zeros(N)
    b[0] = sigma_a * h - c
    b[N - 1] = sigma_a * h - a + 0.5
    b[1:N - 1] = sigma_a * h - c - a

    A = np.zeros(N)
    A[0] = q[0] / b[0]

    B = np.zeros(N)
    B[0] = c / b[0]

    for i in range(N - 1):
        A[i + 1] = (q[i + 1] - a * A[i]) / (b[i + 1] - a * B[i])
        B[i + 1] = c / (b[i + 1] - a * B[i])

    Flow[N - 1] = A[N - 1]
    for i in range(N - 2, -1, -1):
        Flow[i] = A[i] - B[i] * Flow[i + 1]

    k_eff = k_eff * (sigma_gen * np.linalg.norm(Flow)) / (sigma_gen * Last_Norm)

print(k_eff)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
x = np.arange(N)
ax.plot(x, Flow, color='blue', linewidth=1, label="Поток")
plt.xlabel('Расстояние')
plt.ylabel('Поток')
plt.legend()
plt.grid()
plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# x = np.arange(H)
# Flow_Theory = (Q / sigma_a) * (1 - np.cosh(x / L) / np.cosh(H / L))
# ax.plot(x, Flow_Theory, color='blue', linewidth=1, label="Теория")
# ax.plot(x, Flow, color='red', label="Модель")
# plt.xlabel('Шаг пластины')
# plt.ylabel('Поток')
# plt.legend()
# plt.grid()
# plt.show()
