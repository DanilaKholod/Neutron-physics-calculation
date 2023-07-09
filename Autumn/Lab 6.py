import numpy as np
import random
import math
import matplotlib.pyplot as plt

P = 10000
H = 100
S_s = 0.3
S_a = 0.03
S_tot = S_s + S_a
Scattering = 0
Absorption = 0
Dead = 0
Count = 0
S = 0
h = 5
N = int(H / h)

n = np.zeros((3, N))
f = np.zeros((10, N))


def func(x, scattering, absorption, dead, count, s, iter):
    cos = 2 * random.random() - 1
    l = (-1 / S_tot) * np.log(random.random())
    fin = x + cos * l
    if fin >= H:
        dead += 1
        if count != 0:
            s += 1
    else:
        if fin <= 0:
            fin = -fin
        p_s = S_s / S_tot
        p_a = S_a / S_tot
        if random.random() < p_s:
            scattering += 1
            if scattering > 1:
                count += 1
            j = int(fin / h)
            n[0][j] += 1
            f[iter, j] += 1
            scattering, absorption, dead, count, s = func(fin, scattering, absorption, dead, count, s,iter)
        else:
            absorption += 1
            if count != 0:
                s += 1
            j = int(fin / h)
            n[1][j] += 1
    return scattering, absorption, dead, count, s


def avarage(arr, n):
    sum_values = 0.0
    for i in arr:
        sum_values += i
    return sum_values / float(n)


def sigma(arr, n):
    dif = 0.0
    avg = avarage(arr, n)
    for i in arr:
        dif += ((i - avg) ** 2)
    return np.sqrt(dif / (n * (n - 1)))


for j in range(10):
    for i in range(P):
        scattering_i = 0
        absorption_i = 0
        dead_i = 0
        count_i = 0
        s_i = 0
        x = random.random() * H
        (scattering_i, absorption_i, dead_i, count_i, s_i) = func(x, 0, 0, 0, 0, 0, j)
        Count = count_i + Count
        Dead = dead_i + Dead
        Scattering = Scattering + scattering_i
        Absorption = Absorption + absorption_i
        S = S + s_i

print(f"scattering {Scattering} absorption {Absorption} dead {Dead}")
print(n)

for i in range(N):
    n[2, i] = n[0, i] / (S_a * h * P)

avarage_f = np.zeros(N)
sigma_f = np.zeros(N)
for i in range(N):
    avarage_f[i] = np.average(f[:, i])

for i in range(N):
    sigma_f[i] = sigma(f[:, i], 10)


x = np.arange(20)
y = n[2]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, avarage_f, color='teal', linewidth=5)
plt.errorbar(x, avarage_f, yerr=sigma_f)
plt.show()
