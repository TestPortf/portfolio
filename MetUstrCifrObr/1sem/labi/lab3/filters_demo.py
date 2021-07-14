N = 16384
delta_T = 1/44100

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

    # Я дико извиняюсь, но очень надо
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
    import pprint

    filter = []
    for n in range(0, N):
        if n == 0:
            p = 0j
            filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n < N/2:
            p = complex(0, ((2 * np.pi * n)/(N * delta_T)))
            filter.append(transfer_function(p, filter_type, filter_kind, filter_order,  fn, fv))
        elif n == N/2:
            p = abs(complex(0, (np.pi/delta_T)))
            filter.append(transfer_function(p, filter_type, filter_kind, filter_order, fn, fv))
        elif n > N/2:
            filter.append(complex(filter[abs(n - N)].real, -filter[abs(n - N)].imag))

    return filter


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.fftpack as fp

    # Вычисляем частоту
    k = np.linspace(0, int(N/2), int(N/2))
    frequency = k/(N * delta_T)

    filter_type = 'BPF'
    fn = 410
    fv = 390

    # Раскоментите, если захотите получить все порядки всех видов фильтров
    plt.figure()
    for n in range(2, 9):
        our_filter = generate_filter(filter_type, 'Bessel', n, fv=fv, fn=fn)
        plt.plot(frequency, abs(np.array(our_filter))[:len(k)], label='{} порядок фильтра'.format(n))
        plt.title('Фильтры Бесселя')
        plt.legend()

    plt.figure()
    for n in range(2, 9):
        our_filter = generate_filter(filter_type, 'Butterworth', n, fv=fv, fn=fn)
        plt.plot(frequency, abs(np.array(our_filter))[:len(k)], label='{} порядок фильтра'.format(n))
        plt.title('Фильтры Баттерворта')
        plt.legend()
    plt.figure()

    for n in range(2, 9):
        our_filter = generate_filter(filter_type, 'Chebyshev', n, fv=fv, fn=fn)
        plt.plot(frequency, abs(np.array(our_filter))[:len(k)], label='{} порядок фильтра'.format(n))
        plt.title('Фильтры Чебышева')
        plt.legend()

    plt.figure()
    order_number = 3
    plt.title('Сравнение фильтров {} порядка'.format(order_number))
    filter_Bessel = generate_filter(filter_type, 'Bessel', order_number, fv=fv, fn=fn)
    plt.plot(frequency, abs(np.array(filter_Bessel))[:len(k)], label='Фильтр Бесселя'.format(order_number))
    filter_Butterworth = generate_filter(filter_type, 'Butterworth', order_number, fv=fv, fn=fn)
    plt.plot(frequency, abs(np.array(filter_Butterworth))[:len(k)], label='Фильтр Баттерворта'.format(order_number))
    filter_Chebyshev = generate_filter(filter_type, 'Chebyshev', order_number, fv=fv, fn=fn)
    plt.plot(frequency, abs(np.array(filter_Chebyshev))[:len(k)], label='Фильтр Чебышева'.format(order_number))
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
