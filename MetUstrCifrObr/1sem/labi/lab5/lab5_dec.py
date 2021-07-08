# Const
delta_T = 1/1000
N = 65536
q = 8192
from decimal import Decimal


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
            if 'e+' in line:
                int_part, power = line.strip().split('e+')
                array.append(float(int_part)*(10**(int(power))))
            elif 'e-' in line:
                int_part, power = line.strip().split('e-')
                array.append(float(int_part)*(10**(-int(power))))
            else:
                array.append(float(line.strip()))

    return np.array(array)


def signal_plot(input_signal):
    """
    Строит график сигнала и амплитудного спектра

    :input_signal: входной сигнал
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Вычисляем частоту сигнала
    k = np.arange(0, int(q/2), 1)
    frequency = k/(int(q/2) * delta_T)
    
    # Строим сигнал
    plt.figure()
    plt.xlabel('Отсчёты')
    plt.ylabel('Амплитуда')
    plt.plot(input_signal)


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


class ComplexDecimal(object):
    def __init__(self, value):
        self.real = Decimal(value.real)
        self.imag = Decimal(value.imag)

    def __add__(self, other):
        result = ComplexDecimal(self)
        result.real += Decimal(other.real)
        result.imag += Decimal(other.imag)
        return result

    __radd__ = __add__

    def __str__(self):
        return f'({str(self.real)}+{str(self.imag)}j)'

    def sqrt(self):
        result = ComplexDecimal(self)
        if self.imag:
            raise NotImplementedError
        elif self.real > 0:
            result.real = self.real.sqrt()
            return result
        else:
            result.imag = (-self.real).sqrt()
            result.real = Decimal(0)
            return result


def dec_exp(x):
    """
    Return e raised to the power of x.  Result type matches input type.
    """
    from decimal import getcontext

    getcontext().prec += 2
    i, lasts, s, fact, num = 0, 0, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 1
        fact *= i
        num *= x
        s += num / fact
    getcontext().prec -= 2
    return +s


def dec_sinh(x):
    """
    Ищет гиперболический синус. В радианах
    """
    from decimal import getcontext
    
    getcontext().prec += 2
    s = (dec_exp(x) - dec_exp(-x))/2 
    getcontext().prec -= 2
    return +s


def dec_cosh(x):
    """
    Ищет гиперболический синус. В радианах
    """
    from decimal import getcontext
    
    getcontext().prec += 2
    s = (dec_exp(x) + dec_exp(-x))/2 
    getcontext().prec -= 2
    return +s

def dec_pi():
    """
    Compute Pi to the current precision.
    """
    from decimal import getcontext, Decimal

    getcontext().prec += 2  # extra digits for intermediate steps
    three = Decimal(3)      # substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s               # unary plus applies the new precision


def dec_cos(x):
    """
    Return the cosine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).
    """
    from decimal import getcontext, Decimal

    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def dec_sin(x):
    """
    Return the sine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).
    """
    from decimal import getcontext, Decimal

    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def dec_tan(x):
    """
    Возвращает тангенс в радианах
    """
    return dec_sin(x)/dec_cos(x)


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
    from numpy import array
    from decimal import Decimal, getcontext
    
    getcontext().prec = 50 # На всякий случай повышаем количество хранимых данных в переменной

    dec_delta_T = Decimal(delta_T)

    # Вводим константы
    gamma = dec_tan((Decimal(dec_delta_T/2)) * (fv - fn))**(-1)
    sigma = dec_cos((Decimal(dec_delta_T/2)) * (fv + fn))/dec_cos((Decimal(dec_delta_T/2)) * (fv - fn))

    # Определяем вид фильтра 
    if filter_kind == 'Chebyshev':
        epsilon = ((10**Decimal(0.5/10)) - 1).sqrt()
        const_phi = Decimal(1/filter_order) * (Decimal(1/epsilon).ln() + (1 + (1/Decimal(epsilon**2))).sqrt())
        chi = dec_sin(Decimal((2 * k - 1)/(2 * filter_order)) * dec_pi()) * dec_sinh(const_phi)
        mu = dec_sin(Decimal((2 * k - 1)/(2 * filter_order)) * dec_pi())**2 * dec_sinh(const_phi)**2 + \
             dec_cos(Decimal((2 * k - 1)/(2 * filter_order)) * dec_pi())**2 * dec_cosh(const_phi)**2
        if filter_order % 2 == 0 or k <= int(filter_order/2):
            a_constants = [((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)),
                           ((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)) * 2 * chi,
                           ((2**(filter_order - 1)) * epsilon)**(1/int(filter_order/2)) * mu]
        else:
            a_constants = [0, 1, dec_sinh(const_phi)]
    elif filter_kind == 'Butterworth':
        if filter_order % 2 == 0 or k <= int(filter_order/2):
            a_constants = [1, -2 * dec_cos(Decimal((2 * k + filter_order - 1)/(2 * filter_order)) * dec_pi()), 1]
        else:
            a_constants = [0, 1, 1]

    # Определяем тип фильтра
    if filter_type == 'LPF':
        betta = [(fv * dec_delta_T)**2]
        betta.append(2 * betta[0])
        betta.append(betta[0])

        alpha = [a_constants[2] * ((fv * dec_delta_T)**2) + 2 * a_constants[1] * fv * dec_delta_T + 4 * a_constants[0],
                 2 * a_constants[2] * ((fv * dec_delta_T)**2) - 8 * a_constants[0],
                 a_constants[2] * ((fv * dec_delta_T)**2) - 2 * a_constants[1] * fv * dec_delta_T + 4 * a_constants[0]]

        D = (betta[0] * z**2 + betta[1] * z + betta[2])/(alpha[0] * z**2 + alpha[1] * z + alpha[2])

    elif filter_type == 'HPF':
        betta = [4, -8, 4]

        alpha = [a_constants[0] * (fn * dec_delta_T)**2 + 2 * a_constants[1] * fn * dec_delta_T + 4 * a_constants[2],
                 2 * a_constants[0] * (fn * dec_delta_T)**2 - 8 * a_constants[2],
                 a_constants[0] * (fn * dec_delta_T)**2 - 2 * a_constants[1] * fn * dec_delta_T + 4 * a_constants[2]]

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

    return array(D), array(alpha), array(betta)


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
    from decimal import Decimal, getcontext

    getcontext().prec = 50 # На всякий случай повышаем количество хранимых данных в переменной

    # Преобразовываем частоты для правильной работы
    fn = Decimal(fn) * 2 * dec_pi()
    fv = Decimal(fv) * 2 * dec_pi()

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
        alphas.append(current_alpha/current_alpha[0])
        betas.append(current_beta/current_alpha[0])

    return np.array(our_filter), np.array(alphas), np.array(betas)


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
    from scipy.fft import ifft, fft
    from decimal import Decimal, getcontext

    getcontext().prec = 50 # На всякий случай повышаем количество хранимых данных в переменной

    # Импортируем файлы
    #discrete_signal = data_file_to_array('Lab5_Khasanov.dat')
    #bit0 = data_file_to_array('Lab5_Khasanov_s_bit0.dat')
    #bit1 = data_file_to_array('Lab5_Khasanov_s_bit1.dat')
    discrete_signal = data_file_to_array('Lab05_RogozinPl.dat')
    bit0 = data_file_to_array('Lab05_RogozinPl_bit0.dat')
    bit1 = data_file_to_array('Lab05_RogozinPl_bit1.dat')

    # Строим графики сигналов
    signal_plot(abs(fft(discrete_signal)))
    signal_plot(bit0)
    signal_plot(bit1)

    # Вычисляем частоту
    k = np.linspace(0, q, q)
    frequency = k/(q * delta_T)

    # АЧХ сигнала
    plt.figure()
    plt.plot(frequency[:int(q/2)], abs(fft(discrete_signal))[:int(q/2)])
    plt.xscale('log')

    plt.figure()
    plt.plot(frequency[:int(q/2)], abs(fft(discrete_signal))[:int(q/2)])
    plt.xscale('log')

    # Формируем массив коэффициентов
    coefficient = []
    for frequency_number in frequency:
        coefficient.append(complex(0, 2 * np.pi * frequency_number * delta_T)) 

    # Расчёт теоритического рекурсивного фильтра по передаточной функции 
    filter_transfer_function = recursive_transfer_function(np.exp(np.array(coefficient)), 'LPF', 'Chebyshev',
                                                           8, fv=500)

    # Формируем массив коэффициентов высокой точности
    coefficient = []
    for frequency_number in frequency:
        coefficient.append(ComplexDecimal(complex(0, 2 * dec_pi() * Decimal(frequency_number) * Decimal(delta_T))))

    # Генереируем фильтр
    recursive_filter, A, B = generate_recursive_filter(dec_exp(np.array(coefficient)), 'LPF',
                                                             'Butterworth', 4, fv=500)
    print(recursive_filter[1])

    # Строим теоритическое АЧХ фильтра, АЧХ рекурсивного фильтра и АЧХ сигналов, всё в логарифмическом масштабе 
    plt.figure()
    plt.plot(frequency[:int(q/2)], abs(filter_transfer_function)[:int(q/2)], label='теоритически заданный аналоговый прототип')
    plt.plot(frequency[:int(q/2)], abs(recursive_filter)[:int(q/2)], label='рекурсивный фильтр')
    plt.plot(frequency[:int(q/2)], abs(fft(bit0))[:int(q/2)], label='бит0')
    plt.plot(frequency[:int(q/2)], abs(fft(bit1))[:int(q/2)], label='бит1')
    plt.title('Сравнение АЧХ фильтрова')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()

    # Строим теоритическое АЧХ фильтра, АЧХ рекурсивного фильтра и АЧХ сигналов
    plt.figure()
    plt.plot(frequency[:int(q/2)], abs(filter_transfer_function)[:int(q/2)], label='теоритически заданный аналоговый прототип')
    plt.plot(frequency[:int(q/2)], abs(recursive_filter)[:int(q/2)], label='рекурсивный фильтр')
    plt.plot(frequency[:int(q/2)], abs(fft(bit0))[:int(q/2)], label='бит0')
    plt.plot(frequency[:int(q/2)], abs(fft(bit1))[:int(q/2)], label='бит1')
    plt.title('Сравнение АЧХ фильтрова')
    plt.xlabel('Частота, Гц')
    plt.ylabel('Коэффициент подавления')
    plt.legend()

    # Фильтруем сигнал и строим его график
    filtered_signal = signal_filtration_with_recursive_filter(discrete_signal, A, B)
    plt.figure()
    plt.plot(filtered_signal)

    # АЧХ сигнала
    plt.figure()
    plt.plot(abs(fft(filtered_signal))[:int(q/2)])
    plt.xscale('log')

    plt.show()


if __name__ == "__main__":
    main()
