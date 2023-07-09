import numpy as np
import random
import math
import matplotlib.pyplot as plt


# функция решения неоднородного уравнения диффузии
def diffusion(flow, source, diff_len, S_rem, N):
    a = np.zeros(N)
    for i in range(1, N):
        a[i] = - (diff_len[i - 1] + diff_len[i]) / 2

    c = np.zeros(N)
    for i in range(0, N - 1):
        c[i] = - (diff_len[i + 1] + diff_len[i]) / 2

    q = np.zeros(N)
    q[0:N] = source[0:N] * 1

    b = np.zeros(N)
    b[0] = S_rem[0] * 1 - c[0]
    b[N - 1] = S_rem[N - 1] * 1 - a[N - 1] + 0.5
    b[1:N - 1] = S_rem[1:N - 1] * 1 - c[1:N - 1] - a[1:N - 1]

    A = np.zeros(N)
    A[0] = q[0] / b[0]

    B = np.zeros(N)
    B[0] = c[0] / b[0]

    for i in range(N - 1):
        A[i + 1] = (q[i + 1] - a[i + 1] * A[i]) / (b[i + 1] - a[i + 1] * B[i])
        B[i + 1] = c[i + 1] / (b[i + 1] - a[i + 1] * B[i])

    flow[N - 1] = A[N - 1]
    for i in range(N - 2, -1, -1):
        flow[i] = A[i] - B[i] * flow[i + 1]


N = int(120 / 1)
# начальные данные быстрая группа
S_rem_fz = np.zeros(N)
S_rem_fz[0:101] = 1.79E-02  # активная зона
S_rem_fz[101:120] = 5.72E-02  # отражатель

S_gen_fz = np.zeros(N)
S_gen_fz[0:101] = 1.54E-02  # активная зона
S_gen_fz[101:120] = 0  # отражатель

S_s = np.zeros(N)
S_s[0:101] = 1.43E-03  # активная зона
S_s[101:120] = 5.67E-02  # отражатель

D_fz = np.zeros(N)
D_fz[0:101] = 1.07E+00  # активная зона
D_fz[101:120] = 1.33E+00  # отражатель

# начальные данные для тепловой зоны
S_rem_th = np.zeros(N)
S_rem_th[0:101] = 2.49E-01  # активная зона
S_rem_th[101:120] = 1.90E-02  # отражатель

S_gen_th = np.zeros(N)
S_gen_th[0:101] = 4.68E-01  # активная зона
S_gen_th[101:120] = 0  # отражатель

D_th = np.zeros(N)
D_th[0:101] = 5.22E-01  # активная зона
D_th[101:120] = 2.85E-01  # отражатель

# поток быстрой группы нейтронов
Flow_fz = np.zeros(N)
Flow_fz[0:N] = 1

# поток тепловой группы нейтронов
Flow_th = np.zeros(N)
Flow_th[0:N] = 1

Q = np.zeros(N) # создание источника нейтронов
k_eff = 100
k = 1
while not abs(k - k_eff) < 1e-5:
    last_norm = np.linalg.norm(Flow_th * S_gen_th + Flow_fz * S_gen_fz)
    k = k_eff
    Q[0:N] = (S_gen_fz * Flow_fz + S_gen_th * Flow_th) / k_eff
    diffusion(Flow_fz, Q, D_fz, S_rem_fz, N) # решаем неоднородное уравнение диффузии для быстрой группы нейтронов
    diffusion(Flow_th, Flow_fz * S_s, D_th, S_rem_th, N) # решаем неоднородное уравнение диффузии для тепловой группы
    k_eff = k_eff * np.linalg.norm(Flow_th * S_gen_th + Flow_fz * S_gen_fz) / last_norm

print("k эффективное = ", k_eff)

# строим два графика поток
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
x = np.arange(120)
ax.plot(x, Flow_th, color='blue', linewidth=1, label="Тепловая группа")
ax.plot(x, Flow_fz, color='red', label="Быстрая группа")
plt.xlabel('Расстояние')
plt.ylabel('Поток')
plt.legend()
plt.grid()
plt.show()

