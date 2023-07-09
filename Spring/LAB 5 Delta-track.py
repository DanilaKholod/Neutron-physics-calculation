import numpy as np
import math
import random
import matplotlib.pyplot as plt

# задаём количество нейтронов и батчей
number_batch = int(input("Введите количество батчей: "))
number_neutrons = int(input("Введите количество нейтронов: "))
number_neutrons_start = number_neutrons

# задаём параметры сечений
S_gen = 0.03026952
nu_f = 2.4
S_a = 0.03
S_f = S_gen / nu_f
S_s = 0.3  # сечение рассеяния
S_c = S_a - S_gen / nu_f  # сечение захвата
S_tot = S_a + S_s

# задаём количество разбиений
H = 100
h = 5
N = int(H / h)


def func(x, k_eff):
    cos = 2 * random.random() - 1  # направление движения нейтрона
    l = (-1 / S_tot) * np.log(random.random())  # длина пробега
    fin = x + cos * l  # финальная точка для нейтрона
    neutrons_fission = 0
    dead = 0
    if fin >= H:
        dead += 1
    else:
        if fin <= 0:
            fin = -fin
            distation[0] += fin
        p_s = S_s / S_tot  # рассеяние
        p_c = S_c / S_a  # захват
        if random.random() < p_s:  # реакция рассеивания
            position = int(fin / h)
            Number_reactions[0][position] += 1
            if int(x / h) == position:
                distation[position] += abs(fin - x)
            else:
                if fin > x:
                    distation[int(x / h)] += (int(x / h) * h + h) - x
                    for i in range(int(x / h) + 1, position):
                        distation[i] += h
                    distation[position] += fin - position * h
                else:
                    distation[int(x / h)] += x - int(x / h) * h
                    for i in range(position + 1, int(x / h)):
                        distation[i] += h
                    distation[position] += (position * h + h) - fin
            neutrons_fission, fin = func(fin, k_eff)
        else:  # реакия поглощения
            if random.random() < p_c:
                position = int(fin / h)
                Number_reactions[1][position] += 1
            else:  # реакция деления
                p_3 = (nu_f / k_eff)
                neutrons_fission += int(np.floor(p_3))
                p = abs(float(p_3) - np.floor(p_3))
                if random.random() < p:
                    neutrons_fission += 1
    return neutrons_fission, fin


def sigma(array):
    dif = 0.0
    avg = np.average(array)
    for i in array:
        dif += ((i - avg) ** 2)
    return np.sqrt(dif / (len(array) * (len(array) - 1)))


source_x = np.zeros(
    (int(number_batch + 1), int(number_neutrons * 3)))  # кол-во батчей, кол-во нейтронов с запасом в 20% - координаты
# первоначальное распределение нейтронов
for i in range(number_neutrons):
    source_x[0][i] = random.random() * H

number_neutrons_next = 0  # нейтроны следующего поколения
Number_reactions = np.zeros((2, N))  # массив для подсчет актов рассеяния и потока
Flux = np.zeros((number_batch, N))
k_eff = 1
k_eff_last = 1

for batch in range(number_batch):  # цикл для каждого батча
    distation = np.zeros(N)
    for n in range(number_neutrons):  # цикл по количеству нейтронов
        x = source_x[batch][n]  # координата нейтрона - координата деления предыдущего поколения
        neutrons_fission, x_new = func(x, k_eff)
        number_neutrons_next += neutrons_fission
        for i in range(number_neutrons_next - neutrons_fission, number_neutrons_next):
            source_x[batch + 1][i] = x_new
    for position in range(N):
        #Flux[batch][position] = Number_reactions[0, position] / (S_a * h * number_neutrons)
        Flux[batch][position] = distation[position]
    k_eff = k_eff_last * number_neutrons_next / number_neutrons
    k_eff_last = k_eff
    print(k_eff)
    number_neutrons = number_neutrons_next  # следующее поколоение нейтронов - нейтроны деления
    number_neutrons_next = 0  # следующее поколение нейтронов обнуляется

avarage_f = np.zeros(N)
sigma_f = np.zeros(N)
for position in range(N):
    avarage_f[position] = np.average(Flux[int(number_batch * 0.1):, position])

for position in range(N):
    sigma_f[position] = sigma(Flux[int(number_batch * 0.1):, position])

x = np.arange(N)
y = avarage_f

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(x, y, color='teal', linewidth=5)
plt.errorbar(x, y, yerr=sigma_f)
plt.show()
