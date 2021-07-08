# Const
delta_T = 1/1000
N = 1024


def signal_generate(time, harmonics):
    """ 
    Генерирует сигнал с заданными параметрами

    :harmonics: гармоники сигнала в виде словаря {амплитуда:частота} для каждой гармоники
    :time: отсчёты по х-оси
    :returns: сигнал с заданными параметрами
    """
    from numpy import sin, pi, zeros

    generated_signal = zeros(len(time))

    for ind, harmonic in enumerate(harmonics):
        generated_signal += harmonic/(ind + 1) * sin(2 * pi * harmonics[harmonic] * time)

    return generated_signal 


def compare_signals(regular_signal, filtered_signal, plot_range, shift_point):
    """
    Строит график сигнала и его АЧХ. Как бы показывает разницу между сигналами на графике 

    :regular_signal: обычный сигнал
    :filter_signal: фильтрованный сигнал
    :shift_point: Точка сдвига фильтрованного сигнала. Нужна чтобы убрать задержку и было легче сравнить
    сигналы
    plot_range: отрезок на котором будет построен график
    """
    import matplotlib.pyplot as plt
    import numpy as np

    plt.figure()
    plt.plot(plot_range, regular_signal[plot_range], label='Нефильтрованный') 
    plt.plot(plot_range, filtered_signal[np.array(plot_range) + shift_point], label='Фильтрованный')
    plt.xlabel('Отсчёты')
    plt.ylabel('Амплитуда')
    plt.legend()


def autoregressive_spectrum(our_signal, method, order):
    """
    Вычисляет спектр сигнала по авторегрессионной модели

    :our_signal: Сигнал, спектр которого нужно найти
    :method: способ нахождения коэффициентов модели(yule-walker, berg) 
    :order: порядок авторегрессионной модели
    """
    import numpy as np

    # Опеределяем метод
    if method == 'yule-walker':
        # Для начала найдём мат. ожидание сигнала
        expected_value = 0
        for current_signal_value_index in range(0, N):
            expected_value += ((1/N) * our_signal[current_signal_value_index])

        # Найдём ковариационную функцию
        autocovariance_function = []
        for current_autocovariation_number_index in range(0, order + 2):
            sum_value = 0
            for sum_index in range(0, N - order):
                sum_value += ((our_signal[sum_index] - expected_value) * \
                                (our_signal[sum_index + current_autocovariation_number_index] - expected_value))
            autocovariance_function.append(sum_value * 1/N)
            if current_autocovariation_number_index == 342:
                print(2)

        # Вычисляем а-коэффициенты модели по алгоритму Левинсона-Дарбина
        a = np.zeros([order, order])
        a[1][1] = autocovariance_function[1]/autocovariance_function[0]
        for current_a_coefficient_index in range(1, order + 1):
            high_number = 0
            low_number = 0
            for sum_coefficient_index in range(1, order + 1):
                high_number += a[current_a_coefficient_index][sum_coefficient_index] * autocovariance_function[current_a_coefficient_index + 1 - sum_coefficient_index]
                low_number += a[current_a_coefficient_index][sum_coefficient_index] * autocovariance_function[sum_coefficient_index]
            a[current_a_coefficient_index + 1][current_a_coefficient_index + 1] = (autocovariance_function[current_a_coefficient_index + 1] \
                    - high_number)/(autocovariance_function[0] - low_number)

            for coefficient_index in range(1, order + 1):
                a[current_a_coefficient_index + 1][coefficient_index] = a[current_a_coefficient_index][coefficient_index] -\
                        a[current_a_coefficient_index + 1][current_a_coefficient_index + 1] * \
                        a[current_a_coefficient_index][current_a_coefficient_index + 1 - coefficient_index]

        # Вычисляем b коэффициенты
        b = []
        for current_b_coefficient_index in range(1, order + 1):
            sum_value = 0
            for sum_index in range(1, order + 1):
                sum_value += a[current_b_coefficient_index][sum_index] * autocovariance_function[sum_index]
            b.append(np.sqrt(autocovariance_function[0] - sum_value))
        
    # Вычисляем частоту
    k = np.linspace(0, int(N/2), int(N/2))
    frequency = k/(N * delta_T)
    
    # Вычисляем спектральную плотность мощности
    sum_value = 0
    for sum_index in range(1, order + 1):
        # Формируем массив коэффициентов
        coefficient = []
        for frequency_number in frequency:
            coefficient.append(complex(0, 2 * np.pi * frequency_number * delta_T * sum_index))

        sum_value += a[order][sum_index] * np.exp(coefficient)
    spectral_density = delta_T * abs(b[order + 1]/(1 - sum_value))**2    

    return spectral_density



def signal_plot(input_signal):
    """
    Строит график сигнала и амплитудного спектра

    :input_signal: входной сигнал
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Вычисляем частоту сигнала
    k = np.arange(0, int(len(input_signal)/2), 1)
    frequency = k/(int(len(input_signal)) * delta_T)
    
    # Строим сигнал
    plt.figure()
    plt.xlabel('Отсчёты')
    plt.ylabel('Амплитуда')
    plt.plot(frequency, input_signal[:int(N/2)])


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.fft import fft
    import math

    # Задаём характерестики сигнала, который будем генерировать
    A0 = 5
    A1 = 3
    A2 = 4
    
    f0 = 171.36
    f1 = 428.11
    f2 = 362.33

    #signal_harmonics = {A0:f0, A1:f1, A2:f2}
    signal_harmonics = {A0:250}
    signal_discritisation = np.linspace(0, N-1, N) * delta_T

    # Генерируем сигналы
    lambd = 1
    usuall_signal = signal_generate(signal_discritisation, signal_harmonics)
    gaussian_noise_signal = usuall_signal + lambd * np.random.normal(0, 1, N)

    # Вычисляем частоту
    k = np.linspace(0, N, N)
    frequency = k/(N * delta_T)

    signal_plot(gaussian_noise_signal)

    plt.figure()
    plt.plot(frequency[:int(N/2)], autoregressive_spectrum(gaussian_noise_signal, 'yule-walker', int(N/3)))

    plt.show()


if __name__ == "__main__":
    main()
