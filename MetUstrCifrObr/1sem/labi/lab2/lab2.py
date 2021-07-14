# Const
N = 8192
delta_T = 6.832 * 10**(-4)

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

    return array


def calculate_phase_spectrum(complex_spectrum):
    """
    :complex_spectrum: коплексный спектр сигнала
    :returns: фазовый спектр сигнала
    """
    import cmath as cm

    phase_spectrum = []
    for scalar in complex_spectrum:
        phase_spectrum.append(cm.phase(scalar))

    return phase_spectrum

def high_resolution_window_function(step_number):
    """
    :step_number: количество отсчётов
    :returns: окно высокого разрешения
    """
    import numpy as np

    high_resolution_window = []
    for n in range(0, step_number):
        high_resolution_window.append(0.355768 - 0.487396 * np.cos((2 * np.pi * n)/N) + 0.144232 * np.cos((4 * np.pi * n)/N) - 0.012604 *np.cos((6 * np.pi * n)/N))

    return high_resolution_window

def signal_plot(input_signal):
    """
    Строит гарфик сигнала, вычисляет спектры сигнала с заданными оконными функциями

    :input_signal: входной сигнал(дискретный)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.fftpack as fp

    # Вычисляем частоту
    k = np.arange(0, N/2, 1)
    frequency = k/(N * delta_T)
    plt.rcParams['axes.grid'] = True

    # Строим график дискритизированного сигнала
    signal_window = plt.figure()
    signal_ax = signal_window.add_subplot(111)
    signal_ax.plot(range(0, len(input_signal)), input_signal)
    signal_ax.set_title("Дискритизированный сигнал")
    signal_ax.set_xlabel("Отсчёты")
    signal_ax.set_ylabel("Амплитуда, В")

    # Вычисляем амплитудный спектр из БПФ сигнала и строим его график
    fft_sd = fp.fft(input_signal)
    fft_spectrum_window = plt.figure()
    fft_amplitude_spectrum_ax = fft_spectrum_window.add_subplot(211)
    fft_amplitude_spectrum_ax.set_xscale('log')
    fft_amplitude_spectrum_ax.plot(frequency, 20 * np.log10(abs(fft_sd[:len(k)])))
    fft_amplitude_spectrum_ax.set_title("Амплитудный спектр сигнала")
    fft_amplitude_spectrum_ax.set_xlabel("Частота, Гц")
    fft_amplitude_spectrum_ax.set_ylabel("Амплитуда, Дб")
    
    # Вычисляем фазовый спектр из БПФ сигнала и строим его график
    fft_phase_spectrum = calculate_phase_spectrum(fft_sd)
    fft_phase_spectrum_ax = fft_spectrum_window.add_subplot(212)
    fft_phase_spectrum_ax.plot(frequency, fft_phase_spectrum[:len(k)])
    fft_phase_spectrum_ax.set_title("Фазовый спектр сигнала")
    fft_phase_spectrum_ax.set_xlabel("Частота, Гц")
    fft_phase_spectrum_ax.set_ylabel("Фаза, рад")

    # Вычисляем амплитудный спектр сигнала с применением оконной функции высокого разрешения и строим график
    hrw = high_resolution_window_function(N)
    high_res_win_spectrum_window = plt.figure()
    high_res_win_amplitude_spectrum_ax = high_res_win_spectrum_window.add_subplot(211)
    high_res_win_amplitude_spectrum_ax.set_xscale('log')
    high_res_win_amplitude_spectrum_ax.plot(frequency, 20 * np.log10((abs(fp.fft(hrw * input_signal)))[:len(k)]))
    high_res_win_amplitude_spectrum_ax.set_title("Амплитудный спектр сигнала с окном Натталла")
    high_res_win_amplitude_spectrum_ax.set_xlabel("Частота, Гц")
    high_res_win_amplitude_spectrum_ax.set_ylabel("Амплитуда, Дб")

    # Вычисляем фазовый спектр сигнала с применение оконной функции низкого разрешения и строим график
    high_res_win_phase_spectrum = calculate_phase_spectrum(fp.fft(hrw * input_signal))
    high_res_win_phase_spectrum_ax = high_res_win_spectrum_window.add_subplot(212)
    high_res_win_phase_spectrum_ax.plot(frequency, high_res_win_phase_spectrum[:len(k)])
    high_res_win_phase_spectrum_ax.set_title("Фазовый спектр сигнала с окном Натталла")
    high_res_win_phase_spectrum_ax.set_xlabel("Частота, Гц")
    high_res_win_phase_spectrum_ax.set_ylabel("Фаза, рад")

    # Добавляю пространство между графиками
    signal_window.tight_layout()
    fft_spectrum_window.tight_layout()
    high_res_win_spectrum_window.tight_layout()

    signal_window.savefig("Входной сигнал.png")
    fft_spectrum_window.savefig("БПФ спектр сигнала.png")
    high_res_win_spectrum_window.savefig("Оконная функция высокого разрешения.png")


def main():
    import numpy as np
    import matplotlib.pyplot as plt

    # Импортируем сигнал
    file_path = "DSP_Lab_02_File_Var_8_28.dat"
    input_signal = data_file_to_array(file_path)

    # Дополняем сигнал нулями до степени двойки
    while(len(input_signal) != N):
        input_signal.append(0)

    # Строим графики сигнала
    signal_plot(np.array(input_signal))
    plt.show()

if __name__ == "__main__":
    main()
