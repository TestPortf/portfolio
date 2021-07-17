# Const
delta_T = 1/4000
N = 65536
q = 256

def low_resolution_window_function(step_number):
    """
    :step_number: количество отсчётов
    :returns: окно низкого разрешения
    """
    import numpy as np

    low_resolution_window = []
    for n in range(0, step_number):
        low_resolution_window.append(np.exp(-(1/2) * ((n - (q/2))/(0.5 * (q/2)))**2))

    return low_resolution_window


def high_resolution_window_function(step_number):
    """
    :step_number: количество отсчётов
    :returns: окно высокого разрешения
    """
    import numpy as np

    high_resolution_window = []
    for n in range(0, step_number):
        high_resolution_window.append(np.sin((np.pi * n)/q))

    return high_resolution_window


def transfer_function(p, filter_type, filter_kind, filter_order, fn=0, fv=1):
    """
    На основе типа фильтра делает подстановку в передаточную функцию, которую вы написали

    :p: оператор передаточной функции
    :filter_type: Тип фильтра('LPF', 'HPF', 'BPF', 'BSF')
    :filter_kind: Вид фильтра('Chebyshev', 'Butterworth', 'Bessel')
    :filter_order: Порядок фильтра(от 2 до 8)
    :fn: Нижняя граница частоты среза фильтра
    :fv: Верхняя граница частоты среза фильтра
    :returns: расчитанную передаточную функцию

    """
    import numpy as np

    fn = fn * 2 * np.pi
    fv = fv * 2 * np.pi

    # Определяем тип фильтра
    if p == 0:
        changed_p = 0
    else:
        if filter_type == 'LPF':
            changed_p = p/fv
        elif filter_type == 'HPF':
            changed_p = fn/p
        elif filter_type == 'BPF':
            changed_p = ((p**2) + fn * fv)/(p * (fv - fn))
        elif filter_type == 'BSF':
            changed_p = (p * (fv - fn))/((p**2) + fn * fv)
        else:
            print('Неверно указан тип фильтра')

    # Записываем передаточные функции фльтров
    Bessel_transfer_functions = ['3/((changed_p**2) + 3 * changed_p + 3)',
            '15/((changed_p + 2.322) * ((changed_p**2) + 3.678 * changed_p + 6.459))',
            '105/(((changed_p**2) + 5.792 * changed_p + 9.14) * ((changed_p**2) + 4.208 * changed_p + 11.49))',
            '945/((changed_p + 3.647) * ((changed_p**2) + 6.704 * changed_p + 14.27) * ((changed_p**2) + 4.649 * changed_p + 18.16))',
            '10395/(((changed_p**2) + 5.032 * changed_p + 26.51) * ((changed_p**2) + 7.471 * changed_p + 20.85) * ((changed_p**2) + 8.497 * changed_p + 18.80))',
            '135135/((changed_p + 4.972) * ((changed_p**2) + 5.371 * changed_p + 36.6) * ((changed_p**2) + 8.14 * changed_p + 28.94) * ((changed_p**2) + 9.517 * changed_p + 25.67))',
            '2027025/(((changed_p**2) + 5.678 * changed_p + 48.43) * ((changed_p**2) + 8.737 * changed_p + 38.57) * ((changed_p**2) + 10.41 * changed_p + 33.93) * ((changed_p**2) + 11.18 * changed_p + 31.98))']
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
            '0.0447/((changed_p + 0.256) * ((changed_p**2) * 0.114 * changed_p + 1.016) * ((changed_p**2) + 0.319 * changed_p + 0.677) * ((changed_p**2) + 0.462 * changed_p + 0.254))',
            '0.0224/(((changed_p**2) + 0.0872 * changed_p + 1.012) * ((changed_p**2) + 0.248 * changed_p + 0.741) * ((changed_p**2) + 0.372 * changed_p + 0.359) * ((changed_p**2) + 0.439 * changed_p + 0.088))']

    # Определяем вид фильтра и порядок
    if filter_kind == 'Bessel':
        calculated_transfer_function = eval(Bessel_transfer_functions[filter_order - 2])
    elif filter_kind == 'Chebyshev':
        calculated_transfer_function = eval(Chebyshev_transfer_function[filter_order - 2])
    elif filter_kind == 'Butterworth':
        calculated_transfer_function = eval(Butterworth_transfer_function[filter_order - 2])

    return calculated_transfer_function


def generate_filter(filter_type, filter_kind, filter_order, fn=0, fv=1):
    """
    Генерирует фильтр с заданными параметрами

    :filter_type: Тип фильтра('LPF', 'HPF', 'BPF', 'BSF')
    :filter_kind: Вид фильтра('Chebyshev', 'Butterworth', 'Bessel')
    :filter_order: Порядок фильтра(от 2 до 8)
    :fn: Нижняя граница частоты среза фильтра 
    :fv: Верхняя граница частоты среза фильтра
    :returns: Фильтр
    """
    import numpy as np
    import matplotlib.pyplot as plt 

    our_filter = []
    for n in range(0, q):
        if n == 0:
            p = 0j
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n < q/2:
            p = complex(0, ((2 * np.pi * n)/(q * delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order,  fn, fv))
        elif n == q/2:
            p = abs(complex(0, (np.pi/delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n > q/2:
            our_filter.append(complex(our_filter[abs(n - q)].real, -our_filter[abs(n - q)].imag))

    return np.array(our_filter)


def get_nonrecursive_filter(discrete_transfer_function):
    """
    Возвращает нерекурсивный фильтр из передаточной функции

    :discrete_transfer_function: массив дискретной передаточной функции
    :returns: нерекурсивный фильтр
    """
    from scipy.fft import ifft
    import numpy as np
    import matplotlib.pyplot as plt

    fliped_filter_coefficients = ifft(discrete_transfer_function)

    filter_coefficients = np.zeros(q)
    for k in range(0, q):
        if k < q/2:
            filter_coefficients[k] = fliped_filter_coefficients[k + int(q/2)]
        elif k >= q/2:
            filter_coefficients[k] = fliped_filter_coefficients[k - int(q/2)]

    return np.array(filter_coefficients)


def get_nonrecursive_filter_frequency_response(nonrecursive_filter, frequency_for_plot):
    """
    Рисует АЧХ нерекурсивного фильтра

    :nonrecursive_filter: нерекурсивный фильтр
    :frequency_for_plot: частота для подсчёта и построения графика
    """
    import numpy as np
    import matplotlib.pyplot as plt

    freq = range(0, int(delta_T**(-1)))
    frequency_response = []
    for f in freq:
        frequency_response.append(0)
        for k in range(0, q):
            frequency_response[f] += nonrecursive_filter[k] * np.exp(complex(0, -2 * np.pi * f * k * delta_T))
        frequency_response[f] = abs(frequency_response[f])

    return frequency_response


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

def signal_filtration_with_nonrecursive_filter(signal_to_filter, nonrecursive_filter):
    """
    Так как фильтрация с помощью нерекурсивных фильтров это ряд действий. Удобнее вынести этот процесс в
    отдельную функцию.

    :signal_to_filter: сигнал, который нужно отфильтровать
    :nonrecursive_filter: нерекурсивный фильтр
    """
    import numpy as np

    filtered_signal = []
    for n in range(0, N):
        filtered_signal.append(0)
        for k in range(0, q):
            if n - k >= 0:
                filtered_signal[n] += nonrecursive_filter[k] * signal_to_filter[n - k]

    return np.array(filtered_signal)


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.fft import ifft

    # Задаём характерестики сигнала, который будем генерировать
    A0 = 5
    A1 = 3
    A2 = 4
    
    f0 = 171.36
    f1 = 428.11
    f2 = 362.33

    signal_harmonics = {A0:f0, A1:f1, A2:f2}
    signal_discritisation = np.linspace(0, N-1, N) * delta_T

    # Генерируем сигналы
    lambd = 1
    usuall_signal = signal_generate(signal_discritisation, signal_harmonics)
    harmonic_noise_signal = usuall_signal + lambd * signal_generate(signal_discritisation, {1:557, 1:950})
    gaussian_noise_signal = usuall_signal + lambd * np.random.normal(0, 1, N)
   
    # Генереируем фильтр 
    discrete_analog_filter = generate_filter('LPF', 'Chebyshev', 8, fv=500)

    # Вычисляем частоту
    k = np.linspace(0, q, q)
    frequency = k/(q * delta_T)

    # Выводим АЧХ фильтра
    plt.figure()
    plt.plot(frequency[:int(q/2)], abs(discrete_analog_filter)[:int(q/2)])
    plt.title('АЧХ дискритизированного фильтра по передаточной функции')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.savefig('АЧХ дискритизированного фильтра по передаточной функции.png')

    # Генерируем нерекурсивный фильтр и его АЧХ
    nonrecursive_filter = get_nonrecursive_filter(discrete_analog_filter)
    nonrecursive_filter_frequency_response = \
            get_nonrecursive_filter_frequency_response(nonrecursive_filter, frequency)

    # Применяем функцию окна выского разрешения на нерекурсивном фильтре и вычисляем АЧХ полученного фильтра
    nonrecursive_filter_with_window_function = nonrecursive_filter * high_resolution_window_function(q)
    nonrecursive_filter_with_window_function_frequency_response = \
            get_nonrecursive_filter_frequency_response(nonrecursive_filter_with_window_function, frequency)

    # Выводим сравнение нерекурсивных фильтров с применением функции окна и без неё
    plt.figure()
    plt.plot(range(0, int((delta_T**(-1))/2)), 
            nonrecursive_filter_frequency_response[:int((delta_T**(-1))/2)], label='Без функции окна')
    plt.plot(range(0, int((delta_T**(-1))/2)), 
            nonrecursive_filter_with_window_function_frequency_response[:int((delta_T**(-1))/2)], 
            label='С функцией окна')
    plt.title('Сравнение АЧХ фильтров')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.legend()

    # Выводим сравнение нерекурсивных фильтров с применением функции окна и без неё в логарифмических масштабах
    plt.figure()
    plt.plot(range(0, int((delta_T**(-1))/2)), 
            nonrecursive_filter_frequency_response[:int((delta_T**(-1))/2)], label='Без функции окна')
    plt.plot(range(0, int((delta_T**(-1))/2)), 
            nonrecursive_filter_with_window_function_frequency_response[:int((delta_T**(-1))/2)], 
            label='С функцией окна')
    plt.title('Сравнение АЧХ фильтров')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.legend()
    plt.savefig('Сравнение нерекурсивных фильтров с применением функции окна и без.png')
    
    # Выставляем точку смещения фильтрованного сигнала и длинну выборки, которая будет выводиться
    shift_point = 130 #130 для q = 256
    range_const = 200
    plot_range = range(0, range_const)

    # Фильтруем сигнал и строим графики
    filtered_usuall_signal = signal_filtration_with_nonrecursive_filter(usuall_signal, 
            nonrecursive_filter_with_window_function)
    compare_signals(usuall_signal, filtered_usuall_signal, plot_range, shift_point)
    plt.savefig('График фильтрации обычного сигнала.png')

    # Фильтруем сигнал и строим графики
    filtered_gaussian_noise_signal = signal_filtration_with_nonrecursive_filter(gaussian_noise_signal, 
            nonrecursive_filter_with_window_function)
    compare_signals(gaussian_noise_signal, filtered_gaussian_noise_signal, plot_range, shift_point)
    plt.savefig('Результат фильтрации силнала с гауссовским шумом.png')

    # Фильтруем сигнал и строим графики
    filtered_harmonic_noise_signal = signal_filtration_with_nonrecursive_filter(harmonic_noise_signal, 
            nonrecursive_filter_with_window_function)
    compare_signals(harmonic_noise_signal, filtered_harmonic_noise_signal, plot_range, shift_point)
    plt.savefig('Результат фильтрации сигнала с гармоническим шумом.png')

    plt.show()


if __name__ == "__main__":
    main()
