# Const
delta_T = 1/9000
N = 65536
q = 2048


def recursive_transfer_function(z, fn=0, fv=1):
    """
    Генерирует теоритически заданный аналоговый прототип ПФ фильтра Баттерворта
    соответственно входным данным

    :z: оператор передаточной функции
    :fn: Нижняя граница частоты среза фильтра
    :fv: Верхняя граница частоты среза фильтра
    :returns: расчитанную передаточную функцию
    """
    import numpy as np

    fn = fn * 2 * np.pi
    fv = fv * 2 * np.pi

    # Вводим константы
    gamma = 1/np.tan((delta_T/2) * (fv - fn))
    zeta = np.cos((delta_T/2) * (fv + fn))/np.cos((delta_T/2) * (fv - fn))

    # Осуществляем замену переменной р из ФНЧ в ПФ
    changed_p = gamma * (((z**2) - 2 * zeta * z + 1)/((z**2) - 1))

    # Записываем передаточную функцию фильтра Баттерворта 2-ого порядка
    calculated_transfer_function = 1/((changed_p**2) + 1.41421 * changed_p + 1)

    return calculated_transfer_function


def generate_recursive_filter(z, fn=0, fv=1):
    """
    Генерирует рекурсивный фильтр с заданными параметрами и его коэффициенты

    :z: оператор передаточной функции
    :fn: Нижняя граница частоты среза фильтра 
    :fv: Верхняя граница частоты среза фильтра
    :returns: Теоритический фильтр по аналоговой передаточной функции
    :returns: Alpha коэффициенты фильтра
    :returns: Beta коэффициенты фильтра
    """
    import numpy as np

    # Преобразовываем частоты для правильной работы
    fn = fn * 2 * np.pi
    fv = fv * 2 * np.pi

    gamma = 1/np.tan((delta_T/2) * (fv - fn))
    zeta = np.cos((delta_T/2) * (fv + fn))/np.cos((delta_T/2) * (fv - fn))

    # переменные используемые по формуле
    filter_order = 2
    order_number = 1

    # Непосредственно генерируем фильтр
    alphas = []
    betas = []
    our_filter = 1
    a_constants = [1, -2 * np.cos(((2 * order_number + filter_order - 1)/(2 * filter_order)) * np.pi), 1]

    current_betta = [1, 0, -2, 0, 1]

    current_alpha = [gamma**2 * a_constants[0] + gamma * a_constants[1] + a_constants[2],
             -4 * gamma**2 * zeta * a_constants[0] - 2 * gamma * zeta * a_constants[1],
             4 * gamma**2 * zeta**2 * a_constants[0] + 2 * gamma**2 * a_constants[0] - 2 * a_constants[2],
             -4 * gamma**2 * zeta * a_constants[0] + 2 * gamma * zeta * a_constants[1],
             gamma**2 * a_constants[0] - gamma * a_constants[1] + a_constants[2]]

    our_filter = our_filter * ((current_betta[0] * z**4 + current_betta[1] * z**3 + current_betta[2] * z**2 + current_betta[3] * z + current_betta[4])/ \
        (current_alpha[0] * z**4 + current_alpha[1] * z**3 + current_alpha[2] * z**2 + current_alpha[3] * z + current_alpha[4]))
    alphas.append(np.array(current_alpha)/current_alpha[0])
    betas.append(np.array(current_betta)/current_alpha[0])

    return np.array(our_filter), np.array(alphas), np.array(betas)


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


def signal_filtration_with_recursive_filter(signal_to_filter, A_coef, B_coef):
    """
    Фильтрация сигнала с помощью рекурсивного фильтра

    :signal_to_filter: сигнал, который нужно отфильтровать
    :A_coef: Alpha коэффициенты фильтра
    :B_coef: Beta коэффициенты фильтра
    :return: Фильтрованный сигнал
    """
    import numpy as np

    filtered_signal = []
    # Создаём цикл для итерации по каскадам фильтра
    for current_filter_number in range(0, len(A_coef)):
        # Так как фильтр рекурсивный, то нужно организовать использование фильтрованного прошлым каскадом 
        # сигнала начиная со второго по счёту с 1 каскада фильтра
        if current_filter_number != 0:
            signal_to_filter = filtered_signal
            filtered_signal = []
        # Итерация по отсчётам сигнала
        for current_count in range(0, len(signal_to_filter)):
            filtered_signal.append(0)
            # Итариция по коэффициентам beta конкретного каскада
            for current_beta_coefficient in range(0, len(B_coef[current_filter_number])):
                if current_count - current_beta_coefficient >= 0:
                    filtered_signal[current_count] += B_coef[current_filter_number][current_beta_coefficient] * \
                                           signal_to_filter[current_count - current_beta_coefficient]
            # Итерация по коэффициентам alpha конкретного каскада
            for current_alpha_coefficient in range(1, len(A_coef[current_filter_number])):
                if current_count - current_alpha_coefficient >= 0:
                    filtered_signal[current_count] -= A_coef[current_filter_number][current_alpha_coefficient] * \
                                           filtered_signal[current_count - current_alpha_coefficient]

    return np.array(filtered_signal)


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.fft import fft
    import math

    # Задаём характерестики сигнала, который будем генерировать
    A0 = 5
    A1 = 3
    A2 = 4
    
    f0 = 71.36
    f1 = 408.11
    f2 = 400.33

    signal_harmonics = {A1:f1, A2:f2}
    signal_discritisation = np.linspace(0, N-1, N) * delta_T

    # Генерируем сигналы
    lambd = 1
    usuall_signal = signal_generate(signal_discritisation, signal_harmonics)
    harmonic_noise_signal = usuall_signal + lambd * signal_generate(signal_discritisation, {1:557, 1:950})
    gaussian_noise_signal = usuall_signal + lambd * np.random.normal(0, 1, N)

    # Вычисляем частоту
    k = np.linspace(0, q, q)
    frequency = k/(q * delta_T)

    # Формируем массив коэффициентов
    coefficient = []
    for frequency_number in frequency:
        coefficient.append(complex(0, 2 * np.pi * frequency_number * delta_T)) 

    # Расчёт теоритического рекурсивного фильтра по передаточной функции 
    filter_transfer_function = recursive_transfer_function(np.exp(np.array(coefficient)), fv=410, fn=390)

    # Генереируем фильтр
    recursive_filter, A, B = generate_recursive_filter(np.exp(np.array(coefficient)), fv = 410, fn=390)
    print(20*np.log10(abs(generate_recursive_filter(np.exp(complex(0, 2 * np.pi * 39 * delta_T)), fv=410, fn=390)[0])))
    print(20*np.log10(abs(generate_recursive_filter(np.exp(complex(0, 2 * np.pi * 4100 * delta_T)), fv=410, fn=390)[0])))

    plt.figure()
    plt.plot(usuall_signal)

    # Сравниваем теоритическое АЧХ рекурсивного фильтра и АЧХ рассчитаного фльтра 
    plt.figure()
    #plt.plot(frequency, 20 * np.log10(abs(filter_transfer_function)), label='теоритически заданный аналоговый прототип')
    plt.plot(frequency[:int(q/2)], 20 * np.log10(abs(filter_transfer_function))[:int(q/2)], label='теоритически заданный аналоговый прототип')
    plt.plot(frequency[:int(q/2)], 20 * np.log10(abs(recursive_filter))[:int(q/2)], label='рекурсивный фильтр')
    plt.plot(frequency[:int(q/2)], -60 * np.ones(int(q/2)), label='Граница коэффициента подавления -60 дБ')
    plt.title('Сравнение АЧХ фильтрова')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления, дБ')
    plt.xscale('log')
    plt.legend()

    # Выставляем точку смещения фильтрованного сигнала и длинну выборки, которая будет выводиться
    shift_point = 220
    range_start = 0
    range_stop = 1000
    plot_range = range(range_start, range_stop)

    # Фильтруем сигнал без помех и строим графики
    filtered_usuall_signal = signal_filtration_with_recursive_filter(usuall_signal, A, B)
    compare_signals(usuall_signal, filtered_usuall_signal, plot_range, shift_point)

    # Фильтруем сигнал с гармонической помехой и строим графики
    filtered_harmonic_noise_signal = signal_filtration_with_recursive_filter(harmonic_noise_signal, A, B)
    compare_signals(harmonic_noise_signal, filtered_harmonic_noise_signal, plot_range, shift_point)

    # Фильтруем сигнал с гауссовским шумом и строим графики
    filtered_gaussian_noise_signal = signal_filtration_with_recursive_filter(gaussian_noise_signal, A, B)
    compare_signals(gaussian_noise_signal, filtered_gaussian_noise_signal, plot_range, shift_point)
    signal_plot(abs(fft(filtered_gaussian_noise_signal))) 
    plt.xscale('log')

    plt.show()


if __name__ == "__main__":
    main()
