# Const
N = 60000
delta_T = 1/4000
M = 65536

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
    for n in range(0, M):
        if n == 0:
            p = 0j
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n < M/2:
            p = complex(0, ((2 * np.pi * n)/(M * delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order,  fn, fv))
        elif n == M/2:
            p = abs(complex(0, (np.pi/delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n > M/2:
            our_filter.append(complex(our_filter[abs(n - M)].real, -our_filter[abs(n - M)].imag))

    return our_filter


def filter_signal(input_signal):    
    import numpy as np
    import scipy.fft as fp
    
    # Добиваем сигнал до степени двойки
    input_signal = list(input_signal)
    while len(input_signal) != M:
        input_signal.append(0)
    input_signal = np.array(input_signal)

    # Создаём фильтр
    our_filter = generate_filter('LPF', 'Butterworth', 8, fv=500)

    # Применяем его на сигнале
    output_signal = fp.fft(input_signal) * our_filter
    output_signal = fp.ifft(output_signal)

    # Укорачиваем сигнал обратно
    output_signal = list(output_signal)
    while len(output_signal) != N:
        output_signal.pop(len(output_signal) - 1)
    output_signal = np.array(output_signal)

    return output_signal
 
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


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.fft as fp

    A0 = 5
    A1 = 3
    A2 = 4
    
    f0 = 171.36
    f1 = 478.11
    f2 = 362.33

    signal_harmonics = {A0:f0, A1:f1, A2:f2}

    signal_discritisation = np.linspace(0, N-1, N) * delta_T

    # Генерируем сигналы
    lambd = 1
    usuall_signal = signal_generate(signal_discritisation, signal_harmonics)
    harmonic_noise_signal = usuall_signal + lambd * signal_generate(signal_discritisation, {1:557, 1:950})
    gaussian_noise_signal = usuall_signal + lambd * np.random.normal(0, 1, N)
   
    # Вычисляем частоту
    k = np.linspace(0, int(N/2), int(N/2))
    frequency = k/(N * delta_T)

    # Фильтруем сигнал и находим среднеквадратическую ошибку фильтрации
    filtered_usuall_signal = filter_signal(usuall_signal)
    e_usuall_signal = np.sqrt((1/N) * sum((filtered_usuall_signal - usuall_signal)**2))

    # Строим график сигнала
    range_const = 100
    plot_range = range(1, range_const, 1)
    usuall_signal_window = plt.figure()
    usuall_signal_ax = usuall_signal_window.add_subplot(111)
    usuall_signal_ax.plot(plot_range, usuall_signal[plot_range], label='Нефильтрованный')
    usuall_signal_ax.plot(plot_range, filtered_usuall_signal[plot_range], label='Фильтрованный')
    usuall_signal_ax.set_title('Обычный сигнал')
    usuall_signal_ax.set_xlabel('Отсчёты')
    usuall_signal_ax.set_ylabel('Амплитуда')
    usuall_signal_ax.text(0, 0.2, 'Среднеквадратическоая ошибка фильтрации = {}'.format(e_usuall_signal), 
                         fontweight='bold', transform=usuall_signal_ax.transAxes)
    usuall_signal_ax.legend()
    usuall_signal_window.savefig("График фильтрации обычного сигнала.png")

    # Фильтруем сигнал и находим среднеквадратическую ошибку фильтрации
    filtered_gaussian_noise_signal = filter_signal(gaussian_noise_signal)
    e_gaussian_noise_signal = np.sqrt((1/N) * sum((filtered_gaussian_noise_signal - gaussian_noise_signal)**2))

    # Строим график сигнала
    gaussian_noise_siganl_window = plt.figure()
    gaussian_noise_siganl_ax = gaussian_noise_siganl_window.add_subplot(111)
    gaussian_noise_siganl_ax.plot(plot_range, gaussian_noise_signal[plot_range], label='Нефильтрованный') 
    gaussian_noise_siganl_ax.plot(plot_range, filtered_gaussian_noise_signal[plot_range], label='Фильтрованный')
    gaussian_noise_siganl_ax.set_title('Сигнал с белым шумом')
    gaussian_noise_siganl_ax.set_xlabel('Отсчёты')
    gaussian_noise_siganl_ax.set_ylabel('Амплитуда')
    gaussian_noise_siganl_ax.text(0, 0.2, 'Среднеквадратическоая ошибка фильтрации = {}'.format(e_gaussian_noise_signal), 
                         fontweight='bold', transform=gaussian_noise_siganl_ax.transAxes)
    gaussian_noise_siganl_ax.legend()
    gaussian_noise_siganl_window.savefig("Результат фильтрации силнала с гауссовским шумом.png")

    # Фильтруем сигнал и находим среднеквадратическую ошибку фильтрации
    filtered_harmonic_noise_signal = filter_signal(harmonic_noise_signal)
    e_harmonic_noise_signal = np.sqrt((1/N) * sum((harmonic_noise_signal - filtered_harmonic_noise_signal)**2))

    # Строим график сигнала
    harmonic_noise_signal_window = plt.figure()
    harmonic_noise_signal_ax = harmonic_noise_signal_window.add_subplot(111)
    harmonic_noise_signal_ax.plot(plot_range, harmonic_noise_signal[plot_range], label='Нефильтрованный')
    harmonic_noise_signal_ax.plot(plot_range, filtered_harmonic_noise_signal[plot_range], label='Фильтрованный')
    harmonic_noise_signal_ax.set_title('Сигнал с гармоническим шумом')
    harmonic_noise_signal_ax.set_xlabel('Отсчёты')
    harmonic_noise_signal_ax.set_ylabel('Амплитуда')
    harmonic_noise_signal_ax.text(0, 0.2, 'Среднеквадратическоая ошибка фильтрации = {}'.format(e_harmonic_noise_signal), 
                         fontweight='bold', transform=harmonic_noise_signal_ax.transAxes)
    harmonic_noise_signal_ax.legend()
    harmonic_noise_signal_window.savefig("Результат фильтрации сигнала с гармоническим шумом.png")

    # АЧХ фильтра бесселя
    plt.figure()
    plt.plot(frequency, abs(np.array(generate_filter('LPF', 'Bessel', 2, fn=0, fv=500)))[:len(frequency)], label='Бесселя 2-ого порядка')
    plt.plot(frequency, abs(np.array(generate_filter('LPF', 'Butterworth', 2, fn=0, fv=500)))[:len(frequency)], label='Баттерворта')
    plt.plot(frequency, abs(np.array(generate_filter('LPF', 'Chebyshev', 2, fn=0, fv=500)))[:len(frequency)], label='Чебышева')
    plt.plot(frequency, abs(np.array(generate_filter('LPF', 'Bessel', 8, fn=0, fv=500)))[:len(frequency)], label='Бесселя 8-ого порядка')
    plt.title('Сравнение разных видов ФНЧ фильтров')
    plt.legend()
    plt.savefig("АЧХ фильтров")

    plt.show()


if __name__ == "__main__":
    main()
