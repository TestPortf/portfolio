# Const
delta_T = 530 * 10**(-6)
q = 2048 * 4


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
    alphas = np.array(current_alpha)/current_alpha[0]
    betas = np.array(current_betta)/current_alpha[0]

    return np.array(our_filter), np.array(alphas), np.array(betas)


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


def main():
    import numpy as np
    import matplotlib.pyplot as plt

    # Вычисляем частоту
    k = np.linspace(0, q, q - 1)
    frequency = k/(q * delta_T)

    # Формируем массив коэффициентов
    coefficient = []
    for frequency_number in frequency:
        coefficient.append(complex(0, 2 * np.pi * frequency_number * delta_T)) 

    # Расчёт теоритического рекурсивного фильтра по передаточной функции 
    filter_transfer_function = recursive_transfer_function(np.exp(np.array(coefficient)), fv=410, fn=390)

    # Расчёт рекурсивного фильтра
    recursive_filter, A, B = generate_recursive_filter(np.exp(np.array(coefficient)), fv=410, fn=390)
    fn_div_10 = 20 * np.log10(abs(generate_recursive_filter(np.exp(complex(0, 2 * np.pi * 39 * delta_T)), fv=410, fn=390)[0]))
    fv_mul_10 = 20 * np.log10(abs(generate_recursive_filter(np.exp(complex(0, 2 * np.pi * 4100 * delta_T)), fv=410, fn=390)[0]))
    
    # Вывод коэффициентов рекурсивного фильтра
    print("Альфа коэффициенты")
    for coef_number in range(0, len(A)):
        print(f'{coef_number + 1}-ый: {A[coef_number]:.60}')
    print('Бета коэффициенты')
    for coef_number in range(0, len(B)):
        print(f'{coef_number + 1}-ый: {B[coef_number]:.60}')

    # Сравниваем теоритическое АЧХ рекурсивного фильтра и АЧХ рассчитаного фльтра 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(frequency[:int(q/2)], 20 * np.log10(abs(filter_transfer_function))[:int(q/2)], label='AЧХ фильтра рассчитанного по аналоговому прототипу')
    ax.plot(frequency[:int(q/2)], 20 * np.log10(abs(recursive_filter))[:int(q/2)], label='АЧХ фильтра рассчитаного по таблице')
    ax.plot(frequency[:int(q/2)], -60 * np.ones(int(q/2)), '--', label='Граница коэффициента подавления -60 дБ')
    ax.text(0, 0.5, f'Глубина затухания при 39 Гц = {round(fn_div_10)}', transform=ax.transAxes, fontsize=14)
    ax.text(0, 0.4, f'Глубина затухания при 4100 Гц = {round(fv_mul_10)}', transform=ax.transAxes, fontsize=14)
    ax.set_xlabel('Частота, Гц')
    ax.set_ylabel('Коэффициент подавления, дБ')
    ax.set_xscale('log')
    ax.legend()

    #plt.show()


if __name__ == "__main__":
    main()
