import numpy as np
import random
import math
import matplotlib.pyplot as plt

sigma_a = 0.03
sigma_gen = sigma_a
H = 100
h = 1
D = 1
N = int(H / h)
Q = np.zeros(N)

Flow = np.ones(N)

J_input = np.zeros((N, 2))
J_output = np.zeros((N, 2))

Flow_last = np.zeros(N)
k_eff = 100
k = 1
count = 0
while not abs(k - k_eff) < 10e-20 and np.max(np.abs(Flow_last - Flow)) > 10e-20:
    Last_Norm = np.linalg.norm(Flow)
    k = k_eff
    Q[:] = sigma_gen * Flow[:] / k_eff
    Flow_last[:] = Flow[:]
    D = D / h
    A = (6 * D * (1 + 4 * D)) / (1 + 16 * D + 48 * (D ** 2))
    B = (1 - 48 * (D ** 2)) / (1 + 16 * D + 48 * (D ** 2))
    C = (-8 * D) / (1 + 16 * D + 48 * (D ** 2))
    for i in range(N):
        Flow[i] = (Q[i] + (1 - B - C) / h * (J_input[i][1] + J_input[i][0])) / (2 * A / h + sigma_a)
        J_output[i, 1] = B * J_input[i][1] + C * J_input[i][0] + A * Flow[i]
        J_output[i, 0] = C * J_input[i][1] + B * J_input[i][0] + A * Flow[i]
    for i in range(1, N - 1):
        J_input[i][0] = J_output[i - 1][1]
        J_input[i][1] = J_output[i + 1][0]
    J_input[0][0] = J_output[0][0]
    J_input[0][1] = J_output[1][0]
    J_input[N - 1][1] = 0
    J_input[N - 1][0] = J_output[N - 2][1]
    k_eff = k_eff * (sigma_gen * np.linalg.norm(Flow)) / (sigma_gen * Last_Norm)

print(k_eff)

fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
x = np.ones(N)  # координаты середины нода
for i in range(N):
    if i > 0:
        x[i] = x[i - 1] + h
    else:
        x[i] = h / 2

ax1.plot(x, Flow, color='red', linewidth=1)
ax1.legend()

fig1.set_figwidth(12)
fig1.set_figheight(6)
plt.xlabel('Координата')
plt.ylabel('Поток')
plt.grid()
plt.show()
