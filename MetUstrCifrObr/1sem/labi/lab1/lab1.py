def split_signals(signal_amplitudes):
    """ 
    Ищет начала и концы сигналов

    :signal_amplitudes: амплитуды сигнала
    :returns: словарь с началами и концами сигналов

    """
    # Ищем отдельные сигналы
    count = 0
    start_signal = 0
    end_signal = 0
    start_end_siganl_dots = {}
    is_siganl_detected = False
    while True:
        # Ищем начало сигнала
        for index, signal_level in enumerate(signal_amplitudes[end_signal:]):
            if (signal_level > 2) | (signal_level < -2):
                start_signal = index + end_signal - 1
                is_siganl_detected = True
                break
        # Если сигнал не замечен до конца массива, то выходим
        if is_siganl_detected == False:
            break
        # Ищем конец сигнала
        for index, signal_level in enumerate(signal_amplitudes[start_signal:]):
            if (signal_level <= 2) & (signal_level >= -2):
                count += 1
                if count > 70:
                    end_signal = index + start_signal - count + 2
                    count = 0
                    break
            else:
                count = 0
        if is_siganl_detected == True:
            is_siganl_detected = False
        start_end_siganl_dots[start_signal] = end_signal
    return start_end_siganl_dots

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

def main():
    import matplotlib.pyplot as plt
    from statistics import stdev
    import numpy as np

    # Костанты по заданию
    delta_T = 5.888*10**(-5)
    M = 16
    U_max = 6
    file_path = "DSP_Lab_01_File_Var_10_Att_2.dat"

    # Находим ряд расчётных значений
    delta_U = (U_max * 2)/(2**M)
    global_y = data_file_to_array(file_path) 
    global_x = np.array(range(0, len(global_y))) * delta_T
    rms_U = stdev(global_y[:int(len(global_y)/4)] * delta_U)  

    # Построение графика общего сигнала
    global_fig = plt.figure()
    global_ax = global_fig.add_subplot(111)

    global_ax.plot(global_x, global_y * delta_U)

    global_ax.set_xlabel("Время, сек")
    global_ax.set_ylabel("Амплитуда, В")
    global_ax.text(0.3, 0.2, "Среднеквадратическое значение шума АЦП: {}".format(round(rms_U, 7)), fontsize=8, 
                   fontweight="bold", transform=global_ax.transAxes) 
    global_fig.savefig("Общий сигнал.png")

    # Разделение сигнала
    start_end_siganl_dots = split_signals(global_y)

    # Построение графиков сигналов по отдельности с их параметрами
    for start in start_end_siganl_dots:
        tau = start * delta_T
        duration = (start_end_siganl_dots[start] - start) * delta_T
        max_amplitude = max(global_y[start:start_end_siganl_dots[start]]) * delta_U

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(global_x[start:start_end_siganl_dots[start]], global_y[start:start_end_siganl_dots[start]] * delta_U)

        ax.set_xlabel("Время, сек")
        ax.set_ylabel("Амплитуда, В")
        ax.text(0.3, 0.2, "A = {}\n\u03C4 = {}\nT = {}\n".format(max_amplitude, tau, duration), fontsize=14, 
                fontweight="bold", transform=ax.transAxes)
        fig.savefig("Обнаруженная фигура номер {}.png".format(list(start_end_siganl_dots.keys()).index(start)))

    plt.show()

if __name__ == "__main__": 
    main()
