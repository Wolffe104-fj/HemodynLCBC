'''
Hemodynamics of the large circle of blood circulation - HemodynLCBC

by Wolffe104-fj, 2025
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import expit

def heaviside(x):       #Аппроксимация ф. Хевисайда с помощью сигмоиды (для непрерывности)
    return expit(10 * x)

def main():
    #=============== Параметры модели ===============
    #Сопротивления (мм рт.ст. * с / мл)
    R1 = 1.0
    R2 = 0.005
    R3 = 0.013
    R4 = 0.0398

    #Эластичности/емкости (мл / мм рт.ст.)
    C2 = 4.4
    C3 = 1.33
    C4 = 0.8

    L = 0.0005  #Индуктивность (мм рт.ст. * с**2 / мл)
    dt = 0.01   #Шаг интегрирования по времени (с)
    HR = 75     #Частота сердечных сокращений (уд/мин)

    #Пределы изменения эластичности левого желудочка
    Umax = 2.0
    Umin = 0.05

    N_cycles = 10   #Кол-во симулируемых сердечных циклов

    #Начальные значения переменных состояния
    x = np.array([8.0, 7.3, 70.0, 75.0, 20.0])
    #x = [P_желудочек, P_предсердие, P_артерии, P_аорта, Q_поток]

    #=============== Подготовка массивов ===============
    points_per_cycle = int(60 / HR / dt)    #Кол-во временных точек в одном сердечном цикле

    total_points = points_per_cycle * N_cycles  #Общее кол-во временных точек

    #Матрица для хранения значений переменных состояния
    #Каждая строка - одна переменная, каждый столбец - момент времени
    X = np.zeros((len(x), total_points))

    #Массивы для функции упругости и ее производной
    U = np.zeros(total_points)
    dU = np.zeros(total_points)
    h = np.zeros(total_points)

    #Время одного сердечного цикла
    t_cycle = np.linspace(0, 60 / HR - dt, points_per_cycle)

    #=============== Расчет функции упругости ===============
    U[0] = Umin

    for i in range(1, points_per_cycle):
        tn = t_cycle[i] / (0.2 + 0.15 * 60 / HR)    #Нормированное время
        #T = 60/HR - длительность цикла
        #0.2 сек - мин. время активации миокарда (взято среднее)
        #0.15 - доля длительности цикла, добавляемая к базовому времени (взято среднее)

        Un = 1.55 * (tn / 0.7) ** 1.9 / ((1 + (tn / 0.7) ** 1.9) * (1 + (tn / 1.17) ** 21.9))   #Нормированная функция упругости
        #1.55 - нормировочный коэффициент (чтобы max(Un)=1)
        #0.7 - время достижения 50% максимальной активации
        #1.9 - крутизна нарастания активации
        #1.17 - время начала релаксации
        #21.9 - очень большая степень для быстрого спада (релаксация быстрее активации)

        U[i] = (Umax - Umin) * Un + Umin    #Масштабирование функции упругости

        #Производная и вспомогательная функция
        dU[i] = (U[i] - U[i - 1]) / dt
        h[i] = dU[i] / U[i - 1] if U[i - 1] != 0 else 0

    #Копируем функцию упругости для всех циклов (предполагаем периодичность)
    for cycle in range(1, N_cycles):
        start_idx = cycle * points_per_cycle
        end_idx = (cycle + 1) * points_per_cycle
        U[start_idx:end_idx] = U[:points_per_cycle]
        dU[start_idx:end_idx] = dU[:points_per_cycle]
        h[start_idx:end_idx] = h[:points_per_cycle]

    #=============== Интегрирование ===============
    for j in range(N_cycles):
        for i in range(points_per_cycle):

            idx = j * points_per_cycle + i  #Текущий индекс во временном массиве

            x1, x2, x3, x4, x5 = x      #Извлекаем текущие значения переменных

            #Вычисляем условия открытия клапанов через ф. Хевисайда
            H21 = heaviside(x2 - x1)  #Клапан между x2 и x1
            H14 = heaviside(x1 - x4)  #Клапан между x1 и x4

            #Строки соответствуют уравнениям для производных каждой переменной
            A = np.array([
                #Ур-е для x1 (давление в желудочке)
                [h[idx] - U[idx] * H21 / R2 - U[idx] * H14 / R3,
                 U[idx] * H21 / R2,
                 0,
                 U[idx] * H14 / R3,
                 0],

                #Ур-е для x2 (давление в предсердии/у корня аорты)
                [H21 / (R2 * C2),
                 -1 / (R1 * C2) - H21 / (R2 * C2),
                 1 / (R1 * C2),
                 0,
                 0],

                #Ур-е для x3 (давление в артериях)
                [0,
                 1 / (R1 * C3),
                 -1 / (R1 * C3),
                 0,
                 1 / C3],

                #Ур-е для x4 (давление в аорте)
                [H14 / (R3 * C4),
                 0,
                 0,
                 -H14 / (R3 * C4),
                 -1 / C4],

                #Ур-е для x5 (поток через аортальный клапан)
                [0,
                 0,
                 -1 / L,
                 1 / L,
                 -R4 / L]
            ])

            x = dt * A @ x + x      #Интегрирование м. Эйлера

            X[:, idx] = x   #Сохраняем результат

    #=============== Анализ результатов ===============
    #Извлекаем данные последнего сердечного цикла
    last_cycle_start = (N_cycles - 1) * points_per_cycle
    last_cycle_end = N_cycles * points_per_cycle

    #Давление в аорте в последнем цикле
    aortic_pressure = X[3, last_cycle_start:last_cycle_end]

    #Находим систолическое (макс.) давление
    systolic_pressure = np.max(aortic_pressure)

    #Находим диастолическое (мин.) давление
    diastolic_pressure = np.min(aortic_pressure)

    #Пульсовое давление
    pulse_pressure = systolic_pressure - diastolic_pressure

    print("=" * 50)
    print("РЕЗУЛЬТАТЫ МОДЕЛИРОВАНИЯ")
    print("=" * 50)
    print(f"Систолическое давление в аорте: {systolic_pressure:.1f} мм рт.ст.")
    print(f"Диастолическое давление в аорте: {diastolic_pressure:.1f} мм рт.ст.")
    print(f"Пульсовое давление: {pulse_pressure:.1f} мм рт.ст.")
    print(f"Среднее артериальное давление: {np.mean(aortic_pressure):.1f} мм рт.ст.")
    print("=" * 50)

    #=============== Визуализация ===============
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    #Временная ось для последних двух циклов
    time_last_cycles = np.arange(-2 * points_per_cycle, 0) * dt

    #Давление в аорте (последние 2 цикла)
    ax1 = axes[0, 0]
    aortic_last_2cycles = X[3, -2 * points_per_cycle:]
    ax1.plot(time_last_cycles, aortic_last_2cycles, 'b-', linewidth=2)
    ax1.axhline(y=systolic_pressure, color='r', linestyle='--', alpha=0.7,
                label=f'САД: {systolic_pressure:.1f}')
    ax1.axhline(y=diastolic_pressure, color='g', linestyle='--', alpha=0.7,
                label=f'ДАД: {diastolic_pressure:.1f}')
    ax1.set_xlabel('Время (с)')
    ax1.set_ylabel('Давление (мм рт.ст.)')
    ax1.set_title('Давление в аорте (последние 2 цикла)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.fill_between(time_last_cycles, diastolic_pressure, aortic_last_2cycles,
                     alpha=0.2, color='blue')

    #Функция упругости желудочка
    ax2 = axes[0, 1]
    time_full = np.arange(total_points) * dt
    ax2.plot(time_full[:points_per_cycle * 2], U[:points_per_cycle * 2], 'r-', linewidth=2)
    ax2.set_xlabel('Время (с)')
    ax2.set_ylabel('Упругость')
    ax2.set_title('Функция упругости левого желудочка')
    ax2.grid(True, alpha=0.3)
    ax2.fill_between(time_full[:points_per_cycle * 2], Umin, U[:points_per_cycle * 2],
                     alpha=0.2, color='red')

    #Все давления в последнем цикле
    ax3 = axes[1, 0]
    time_last_cycle = np.arange(points_per_cycle) * dt
    labels = ['Желудочек', 'Предсердие', 'Артерии', 'Аорта']
    colors = ['red', 'orange', 'green', 'blue']

    for idx in range(4):
        ax3.plot(time_last_cycle, X[idx, last_cycle_start:last_cycle_end],
                 label=labels[idx], color=colors[idx], linewidth=2)

    ax3.set_xlabel('Время (с)')
    ax3.set_ylabel('Давление (мм рт.ст.)')
    ax3.set_title('Давления в различных отделах (последний цикл)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    #Поток крови
    ax4 = axes[1, 1]
    flow = X[4, last_cycle_start:last_cycle_end]
    ax4.plot(time_last_cycle, flow, 'purple', linewidth=2)
    ax4.set_xlabel('Время (с)')
    ax4.set_ylabel('Поток (мл/с)')
    ax4.set_title('Поток крови через аортальный клапан')
    ax4.grid(True, alpha=0.3)
    ax4.fill_between(time_last_cycle, 0, flow, alpha=0.2, color='purple')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
