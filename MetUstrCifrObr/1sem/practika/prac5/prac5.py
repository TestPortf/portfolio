# Const
delta_T = 1/4000
N = 65536
q = 256


def recursive_transfer_function(z, filter_type, filter_kind, filter_order, fn=0, fv=1):
    """
    Генерирует теоритически заданный аналоговый прототип фильтра соответственно входным данным

    :z: оператор передаточной функции
    :filter_type: Тип фильтра('LPF', 'HPF', 'BPF', 'BSF')
    :filter_kind: Вид фильтра('Chebyshev', 'Butterworth')
    :filter_order: Порядок фильтра(от 2 до 8)
    :fn: Нижняя граница частоты среза фильтра
    :fv: Верхняя граница частоты среза фильтра
    :returns: расчитанную передаточную функцию
    """
    import numpy as np

    fn = fn * 2 * np.pi
    fv = fv * 2 * np.pi

    # Вводим константы
    gamma = 1/np.tan((delta_T/2) * (fv - fn))
    sigma = np.cos((delta_T/2) * (fv + fn))/np.cos((delta_T/2) * (fv - fn))

    # Определяем тип фильтра
    if filter_type == 'LPF':
        changed_p = (2/(fv * delta_T)) * ((z - 1)/(z + 1))
    elif filter_type == 'HPF':
        changed_p = ((fn * delta_T)/2) * ((z + 1)/(z - 1))
    elif filter_type == 'BPF':
        changed_p = gamma * (((z**2) - 2 * sigma * z + 1)/((z**2) - 1))
    elif filter_type == 'BSF':
        changed_p = ((z**2) - 1)/(gamma * (((z**2) - 2 * sigma * z + 1)))
    else:
        print('Неверно указан тип фильтра')

    # Записываем передаточные функции фльтров
    Butterworth_transfer_function = ['1/((changed_p**2) + 1.41421 * changed_p + 1)',
            '1/((changed_p + 1) * ((changed_p**2) + changed_p + 1))',
            '1/(((changed_p**2) + 1.84776 * changed_p + 1) * ((changed_p**2) + 0.76537 * changed_p + 1))',
            '1/((changed_p + 1) * ((changed_p**2) + 1.61803 * changed_p + 1) * ((changed_p**2) + 0.61803 * changed_p + 1))',
            '1/(((changed_p**2) + 1.93185 * changed_p + 1) * ((changed_p**2) + 1.41421 * changed_p + 1) * ((changed_p**2) + 0.51764 * changed_p + 1))',
            '1/((changed_p + 1) * ((changed_p**2) + 1.80194 * changed_p + 1) * ((changed_p**2) + 1.24698 * changed_p + 1) * ((changed_p**2) + 0.44504 * changed_p + 1))',
            '1/(((changed_p**2) + 1.96157 * changed_p + 1) * ((changed_p**2) + 1.66294 * changed_p + 1) * ((changed_p**2) + 1.11114 * changed_p + 1) * ((changed_p**2) + 0.39018 * changed_p + 1))']
    Chebyshev_transfer_function = ['1.431/((changed_p**2) + 1.426 * changed_p + 1.516)',
            '0.716/((changed_p + 0.626) * ((changed_p**2) + 0.626 * changed_p + 1.142))',
            '0.358/(((changed_p**2) + 0.351 * changed_p + 1.064) * ((changed_p**2) + 0.847 * changed_p + 0.356))',
            '0.1789/((changed_p + 0.362) * ((changed_p**2) + 0.224 * changed_p + 1.036) * ((changed_p**2) + 0.586 * changed_p + 0.477))',
            '0.0895/(((changed_p**2) + 0.155 * changed_p + 1.023) * ((changed_p**2) + 0.424 * changed_p + 0.59) * ((changed_p**2) + 0.58 * changed_p + 0.157))',
            '0.0447/((changed_p + 0.256) * ((changed_p**2) + 0.114 * changed_p + 1.016) * ((changed_p**2) + 0.319 * changed_p + 0.677) * ((changed_p**2) + 0.462 * changed_p + 0.254))',
            '0.0224/(((changed_p**2) + 0.0872 * changed_p + 1.012) * ((changed_p**2) + 0.248 * changed_p + 0.741) * ((changed_p**2) + 0.372 * changed_p + 0.359) * ((changed_p**2) + 0.439 * changed_p + 0.088))']

    # Определяем вид фильтра и порядок
    if filter_kind == 'Chebyshev':
        calculated_transfer_function = eval(Chebyshev_transfer_function[filter_order - 2])
    elif filter_kind == 'Butterworth':
        calculated_transfer_function = eval(Butterworth_transfer_function[filter_order - 2])

    return calculated_transfer_function


def generate_filter_chain(z, k, filter_type, filter_kind, filter_order, fn=0, fv=1):
    """
    Генерирует звено рекурсивного филтра с заданными параметрами и коэффициенты звена

    :z: оператор передаточной функции
    :k: порядок звена
    :filter_type: Тип фильтра('LPF', 'HPF', 'BPF', 'BSF')
    :filter_kind: Вид фильтра('Chebyshev', 'Butterworth')
    :filter_order: Порядок фильтра(от 2 до 8)
    :fn: Нижняя граница частоты среза фильтра
    :fv: Верхняя граница частоты среза фильтра
    :returns: расчитанную формулу сомножителей фильтра
    :returns: Alpha коэффициент звена
    :returns: Beta коэффициент звена
    """
    import numpy as np

    # Вводим константы
    gamma = 1/np.tan((delta_T/2) * (fv - fn))
    sigma = np.cos((delta_T/2) * (fv + fn))/np.cos((delta_T/2) * (fv - fn))

    # Определяем вид фильтра 
    if filter_kind == 'Chebyshev':
        epsilon = np.sqrt((10**(0.5/10)) - 1)
        const_phi = (1/filter_order) * np.log((1/epsilon) + np.sqrt(1 + (1/(epsilon**2))))
        chi = np.sin((2 * k - 1)/(2 * filter_order) * np.pi) * np.sinh(const_phi)
        mu = np.sin((2 * k - 1)/(2 * filter_order) * np.pi)**2 * np.sinh(const_phi)**2 + \
             np.cos((2 * k - 1)/(2 * filter_order) * np.pi)**2 * np.cosh(const_phi)**2
        if filter_order % 2 == 0 or k <= int(filter_order/2):
            a_constants = [((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)), ((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)) * 2 * chi,
                           ((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)) * mu]
        else:
            a_constants = [0, 1, np.sinh(const_phi)]
    elif filter_kind == 'Butterworth':
        if filter_order % 2 == 0 or k <= int(filter_order/2):
            a_constants = [1, -2 * np.cos(((2 * k + filter_order - 1)/(2 * filter_order)) * np.pi), 1]
        else:
            a_constants = [0, 1, 1]

    # Определяем тип фильтра
    if filter_type == 'LPF':
        betta = [(fv * delta_T) ** 2]
        betta.append(2 * betta[0])
        betta.append(betta[0])

        alpha = [a_constants[2] * ((fv * delta_T) ** 2) + 2 * a_constants[1] * fv * delta_T + 4 * a_constants[0],
                 2 * a_constants[2] * ((fv * delta_T) ** 2) - 8 * a_constants[0],
                 a_constants[2] * ((fv * delta_T) ** 2) - 2 * a_constants[1] * fv * delta_T + 4 * a_constants[0]]

        D = (betta[0] * z**2 + betta[1] * z + betta[2])/(alpha[0] * z**2 + alpha[1] * z + alpha[2])

    elif filter_type == 'HPF':
        betta = [4, -8, 4]

        alpha = [a_constants[0] * (fn * delta_T)**2 + 2 * a_constants[1] * fn * delta_T + 4 * a_constants[2],
                 2 * a_constants[0] * (fn * delta_T)**2 - 8 * a_constants[2],
                 a_constants[0] * (fn * delta_T)**2 - 2 * a_constants[1] * fn * delta_T + 4 * a_constants[2]]

        D = (betta[0] * z**2 + betta[1] * z + betta[2])/(alpha[0] * z**2 + alpha[1] * z + alpha[2])

    elif filter_type == 'BPF':
        betta = [1, 0, -2, 0, 1]

        alpha = [gamma**2 * a_constants[0] + gamma * a_constants[1] + a_constants[2],
                 -4 * gamma**2 * sigma * a_constants[0] - 2 * gamma * sigma * a_constants[1],
                 4 * gamma**2 * sigma**2 * a_constants[0] + 2 * gamma**2 * a_constants[0] - 2 * a_constants[2],
                 -4 * gamma**2 * sigma * a_constants[0] + 2 * gamma * sigma * a_constants[1],
                 gamma**2 * a_constants[0] - gamma * a_constants[1] + a_constants[2]]

        D = (betta[0] * z**4 + betta[1] * z**3 + betta[2] * z**2 + betta[3] * z + betta[4])/ \
            (alpha[0] * z**4 + alpha[1] * z**3 + alpha[2] * z**2 + alpha[3] * z + alpha[4])

    elif filter_type == 'BSF':
        betta = [gamma**2, -4 * gamma**2 * sigma, 4 * gamma**2 * sigma**2 + 2 * gamma**2]
        betta.append(betta[1])
        betta.append(betta[0])

        alpha = [gamma**2 * a_constants[2] + gamma * a_constants[1] + a_constants[0],
                 -4 * gamma**2 * sigma * a_constants[2] - 2 * gamma * sigma * a_constants[1],
                 4 * gamma**2 * sigma**2 * a_constants[2] + 2 * gamma**2 * a_constants[2] - 2 * a_constants[0],
                 -4 * gamma**2 * sigma * a_constants[2] + 2 * gamma * sigma * a_constants[1],
                 gamma**2 * a_constants[2] - gamma * a_constants[1] + a_constants[0]]

        D = (betta[0] * z**4 + betta[1] * z**3 + betta[2] * z**2 + betta[3] * z + betta[4])/ \
            (alpha[0] * z**4 + alpha[1] * z**3 + alpha[2] * z**2 + alpha[3] * z + alpha[4])

    else:
        print('Неверно указан тип фильтра')

    return np.array(D), alpha, betta


def generate_recursive_filter(z, filter_type, filter_kind, filter_order, fn=0, fv=1):
    """
    Генерирует рекурсивный фильтр с заданными параметрами и его коэффициенты

    :z: оператор передаточной функции
    :filter_type: Тип фильтра('LPF', 'HPF', 'BPF', 'BSF')
    :filter_kind: Вид фильтра('Chebyshev', 'Butterworth')
    :filter_order: Порядок фильтра(от 2 до 8)
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

    # Обрабатываем нечётность фильтра
    order_modificator = 0
    if filter_order % 2 == 1:
        order_modificator = 1

    # Непосредственно генерируем фильтр
    alphas = []
    betas = []
    our_filter = 1
    for order_number in range(1, int((filter_order + order_modificator)/2) + 1):
        filter_coef, current_alpha, current_beta = generate_filter_chain(z, order_number, filter_type, filter_kind,
                                                                         filter_order, fn=fn, fv=fv)
        our_filter = our_filter * filter_coef
        alphas.append(np.array(current_alpha)/current_alpha[0])
        betas.append(np.array(current_beta)/current_alpha[0])

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
    print(len(input_signal))
    print(len(frequency))
    
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
    
    f0 = 171.36
    f1 = 428.11
    f2 = 362.33

    signal_harmonics = {A0:f0, A1:f1, A2:f2}
    #signal_harmonics = {A0:250}
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
    filter_transfer_function = recursive_transfer_function(np.exp(np.array(coefficient)), 'LPF', 'Butterworth',
                                                           8, fv=500)

    # Генереируем фильтр
    recursive_filter, A, B = generate_recursive_filter(np.exp(np.array(coefficient)), 'LPF',
                                                             'Butterworth', 8, fv = 500)

    plt.figure()
    plt.plot(usuall_signal)

    # Сравниваем теоритическое АЧХ рекурсивного фильтра и АЧХ рассчитаного фльтра 
    plt.figure()
    plt.plot(frequency[:int(q/2)], 20 * np.log10(abs(filter_transfer_function))[:int(q/2)], label='теоритически заданный аналоговый прототип')
    plt.plot(frequency[:int(q/2)], 20 * np.log10(abs(recursive_filter))[:int(q/2)], label='рекурсивный фильтр')
    plt.title('Сравнение АЧХ фильтров')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.xscale('log')
    plt.legend()
    plt.savefig('Сравнение АЧХ фильтров.png')

    # Выставляем точку смещения фильтрованного сигнала и длинну выборки, которая будет выводиться
    shift_point = 0
    range_start = 0
    range_stop = 200
    plot_range = range(range_start, range_stop)

    # Фильтруем сигнал без помех и строим графики
    filtered_usuall_signal = signal_filtration_with_recursive_filter(usuall_signal, A, B)
    compare_signals(usuall_signal, filtered_usuall_signal, plot_range, shift_point)
    plt.savefig('График фильтрации обычного сигнала.png')

    # Фильтруем сигнал с гармонической помехой и строим графики
    filtered_harmonic_noise_signal = signal_filtration_with_recursive_filter(harmonic_noise_signal, A, B)
    compare_signals(harmonic_noise_signal, filtered_harmonic_noise_signal, plot_range, shift_point)
    plt.savefig('Результат фильтрации сигнала с гармоническим шумом.png')

    # Фильтруем сигнал с гауссовским шумом и строим графики
    filtered_gaussian_noise_signal = signal_filtration_with_recursive_filter(gaussian_noise_signal, A, B)
    compare_signals(gaussian_noise_signal, filtered_gaussian_noise_signal, plot_range, shift_point)
    plt.savefig('Результат фильтрации силнала с гауссовским шумом.png')

    signal_plot(abs(fft(filtered_gaussian_noise_signal))) 
    plt.xscale('log')
    plt.savefig('Спектр сигнала после фильтрации.png')

    plt.show()


if __name__ == "__main__":
    main()
