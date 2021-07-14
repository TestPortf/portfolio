# Const
N = 14947
delta_T = 1/44100

def data_file_to_array(file_name):
    """
    Преобразует данные из файла в массив амплитуд

    :file_name: путь к файлу в формате строки
    :returns: массив амплитуд

    """
    import numpy as np

    array = []
    with open(file_name, "r") as file:
        for line in file:
            if 'e' in line:
                int_part, power = line.strip().split('e+')
                array.append(float(int_part)*(10**(int(power))))
            else:
                array.append(float(line.strip()))

    return np.array(array)

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
            '0.0447/((changed_p + 0.256) * ((changed_p**2) + 0.114 * changed_p + 1.016) * ((changed_p**2) + 0.319 * changed_p + 0.677) * ((changed_p**2) + 0.462 * changed_p + 0.254))',
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
    for n in range(0, N):
        if n == 0:
            p = 0j
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n < N/2:
            p = complex(0, ((2 * np.pi * n)/(N * delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order,  fn, fv))
        elif n == N/2:
            p = abs(complex(0, (np.pi/delta_T)))
            our_filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n > N/2:
            our_filter.append(complex(our_filter[abs(n - N)].real, -our_filter[abs(n - N)].imag))

    return our_filter


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.fftpack as fp

    # Импортируем данные из файла в массив
    file_name = 'Lab_3_28.dat'
    usuall_signal = data_file_to_array(file_name)

    # Вычисляем частоту
    k = np.linspace(0, int(N/2), int(N/2))
    frequency = k/(N * delta_T)

    # Строим график сигнала и его АЧХ
    signal_window = plt.figure()
    signal_ax = signal_window.add_subplot(211)
    signal_ax.plot(range(len(usuall_signal)), usuall_signal)
    amplitude_spectum_ax = signal_window.add_subplot(212)
    amplitude_spectum_ax.plot(frequency, abs(fp.fft(usuall_signal))[:len(k)])
    signal_window.savefig("Сигнал и его АЧХ.png")

    # Задаём и применяем фильтры
    first_filter = generate_filter('BSF', 'Bessel', 2, fv=43, fn=69)
    second_filter = generate_filter('BSF', 'Bessel', 2,  fv=383.5, fn=403.5)
    third_filter = generate_filter('BSF', 'Bessel', 2, fv=3860, fn=3920)
    filtered_signal = fp.fft(usuall_signal) * first_filter * second_filter * third_filter

    # Строим график фильтрованного сигнала и его АЧХ
    filtered_signal_window = plt.figure()
    signal_ax = filtered_signal_window.add_subplot(211)
    signal_ax.plot(range(len(filtered_signal)), fp.ifft(filtered_signal))
    amplitude_spectum_ax = filtered_signal_window.add_subplot(212)
    amplitude_spectum_ax.plot(frequency, abs(filtered_signal)[:len(k)])
    filtered_signal_window.savefig("Фильтрованный сигнал и его АЧХ.png")

    plt.show()


if __name__ == "__main__":
    main()
