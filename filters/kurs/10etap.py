# Const
delta_T = 1/9000
N = 2000
q = 3


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
        generated_signal += harmonic * sin(2 * pi * harmonics[harmonic] * time)

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
    # Итерация по отсчётам сигнала
    for current_count in range(0, len(signal_to_filter)):
        filtered_signal.append(0)
        # Итариция по коэффициентам beta конкретного каскада
        for current_beta_coefficient in range(0, len(B_coef)):
            if current_count - current_beta_coefficient >= 0:
                filtered_signal[current_count] += B_coef[current_beta_coefficient] * \
                                       signal_to_filter[current_count - current_beta_coefficient]
        # Итерация по коэффициентам alpha конкретного каскада
        for current_alpha_coefficient in range(1, len(A_coef)):
            if current_count - current_alpha_coefficient >= 0:
                filtered_signal[current_count] -= A_coef[current_alpha_coefficient] * \
                                       filtered_signal[current_count - current_alpha_coefficient]

    return np.array(filtered_signal)

def main():
    import numpy as np
    import matplotlib.pyplot as plt

    # Задаём характерестики сигнала, который будем генерировать
    A1 = 3
    A2 = 2

    f1 = 200
    f2 = 400.33

    signal_discritisation = np.linspace(0, N-1, N) * delta_T

    # Генерируем сигналы
    lambd = 1
    signal_in_bandwith = signal_generate(signal_discritisation, {A2:f2})
    signal_out_of_bandwith = signal_generate(signal_discritisation, {A1:f1})
    gaussian_noise_signal = signal_in_bandwith + lambd * np.random.normal(0, 1, N)

    # Вычисляем частоту
    k = np.linspace(0, q, q)
    frequency = k/(q * delta_T)

    # Генереируем фильтр
    A_coef = np.array([1, -3.8261594180878564230852134642191231250762939453125, 5.64021677662502707306657612207345664501190185546875, 
-3.788568424693382841184075005003251135349273681640625, 
0.98044753188059352577710114928777329623699188232421875])
    B_coef = np.array([0.0000482615212977725109556002835997645661336719058454036712646484375, 0, -0.000096523042595545021911200567199529132267343811690807342529296875,0, 0.0000482615212977725109556002835997645661336719058454036712646484375])
    
    xk = 0
    xk_1 = 0
    xk_2 = 0
    xk_3 = 0
    xk_4 = 0

    yk = 0
    yk_1 = 0
    yk_2 = 0
    yk_3 = 0
    yk_4 = 0

    filtered_gaussian_noise_signal = []
    for number in range(0, len(gaussian_noise_signal)):            
        xk_4 = xk_3
        xk_3 = xk_2
        xk_2 = xk_1
        xk_1 = xk

        xk = gaussian_noise_signal[number]

        yk_4 = yk_3
        yk_3 = yk_2
        yk_2 = yk_1
        yk_1 = yk
        yk = B_coef[0] * xk + B_coef[2] * xk_2 + B_coef[4] * xk_4 - A_coef[1] * yk_1 - A_coef[2] * yk_2 - A_coef[3] * yk_3 - A_coef[4] * yk_4
        filtered_gaussian_noise_signal.append(yk)
    filtered_gaussian_noise_signal = np.array(filtered_gaussian_noise_signal)


    # Выставляем точку смещения фильтрованного сигнала и длинну выборки, которая будет выводиться
    shift_point = 220
    range_start = 0
    range_stop = 1000
    plot_range = range(range_start, range_stop)

    # # Фильтруем сигнал без помех и строим графики
    # filtered_usuall_signal = signal_filtration_with_recursive_filter(signal_in_bandwith, A, B)
    # compare_signals(signal_in_bandwith, filtered_usuall_signal, plot_range, shift_point)

    # # Фильтруем сигнал с гармонической помехой и строим графики
    # filtered_harmonic_noise_signal = signal_filtration_with_recursive_filter(signal_out_of_bandwith, A, B)
    # compare_signals(signal_out_of_bandwith, filtered_harmonic_noise_signal, plot_range, shift_point)

    # Фильтруем сигнал с гауссовским шумом и строим графики
    #filtered_gaussian_noise_signal = signal_filtration_with_recursive_filter(gaussian_noise_signal, A, B)
    compare_signals(gaussian_noise_signal, filtered_gaussian_noise_signal, plot_range, shift_point)

    filtered_gaussian_noise_signal = signal_filtration_with_recursive_filter(gaussian_noise_signal, A_coef, B_coef)
    compare_signals(gaussian_noise_signal, filtered_gaussian_noise_signal, plot_range, shift_point)

    plt.show()


if __name__ == "__main__":
    main()
