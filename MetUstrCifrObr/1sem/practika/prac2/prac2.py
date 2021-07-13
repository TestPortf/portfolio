# Const
N = 1024
delta_T = 1/4000

A0 = 0.01
A1 = 5
A2 = 3

f0 = 662.33
f1 = 171.36
f2 = 173.11

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

def low_resolution_window_function(step_number):
    """
    :step_number: количество отсчётов
    :returns: окно низкого разрешения 
    """
    import numpy as np

    low_resolution_window = []
    for n in range(0, step_number):
        low_resolution_window.append(np.exp(-(1/2) * ((n - (N/2))/(0.5 * (N/2)))**2))
    
    return low_resolution_window

def high_resolution_window_function(step_number):
    """
    :step_number: количество отсчётов
    :returns: окно высокого разрешения
    """
    import numpy as np

    high_resolution_window = []
    for n in range(0, step_number):
        high_resolution_window.append(np.sin((np.pi * n)/N))

    return high_resolution_window

def signal_plot(input_signal, signal_name):
    """
    Строит гарфик сигнала, вычисляет спектры сигнала с заданными оконными функциями

    :input_signal: входной сигнал(дискретный)
    :signal_name: название сигнала, которое будет записано в файл
    """
    import numpy as np
    import cmath as cm
    import matplotlib.pyplot as plt
    import scipy.fftpack as fp
    import scipy.interpolate as ip

    # Вычисляем частоту
    k = np.linspace(0, int(N/2), int(N/2))
    frequency = k/(N * delta_T)

    # Строим график дискретизированного сигнала
    signal_window = plt.figure()
    signal_ax = signal_window.add_subplot(111)
    signal_ax.plot(range(0, len(input_signal)), input_signal)
    signal_ax.set_title("Дискретизированный сигнал")
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

    # Вычисляем амплитудный спектр сигнала с применением оконной функции низкого разрешения и строим график
    lrw = low_resolution_window_function(N)
    low_res_win_spectrum_window = plt.figure()
    low_res_win_amplitude_spectrum_ax = low_res_win_spectrum_window.add_subplot(211)
    low_res_win_amplitude_spectrum_ax.set_xscale('log')
    low_res_win_amplitude_spectrum_ax.plot(frequency, 20 * np.log10((abs(fp.fft(lrw * input_signal)))[:len(k)]))
    low_res_win_amplitude_spectrum_ax.set_title("Амплитудный спектр сигнала с окном низкого разрешения")
    low_res_win_amplitude_spectrum_ax.set_xlabel("Частота, Гц")
    low_res_win_amplitude_spectrum_ax.set_ylabel("Амплитуда, Дб")

    # Вычисляем фазовый спектр сигнала с применение оконной функции низкого разрешения и строим график
    low_res_win_phase_spectrum = calculate_phase_spectrum(fp.fft(lrw * input_signal))
    low_res_win_phase_spectrum_ax = low_res_win_spectrum_window.add_subplot(212)
    low_res_win_phase_spectrum_ax.plot(frequency, low_res_win_phase_spectrum[:len(k)])
    low_res_win_phase_spectrum_ax.set_title("Фазовый спектр сигнала с окном низкого разрешения")
    low_res_win_phase_spectrum_ax.set_xlabel("Частота, Гц")
    low_res_win_phase_spectrum_ax.set_ylabel("Фаза, рад")

    # Вычисляем амплитудный спектр сигнала с применением оконной функции высокого разрешения и строим график
    hrw = high_resolution_window_function(N)
    high_res_win_spectrum_window = plt.figure()
    high_res_win_amplitude_spectrum_ax = high_res_win_spectrum_window.add_subplot(211)
    high_res_win_amplitude_spectrum_ax.set_xscale('log')
    high_res_win_amplitude_spectrum_ax.plot(frequency, 20 * np.log10((abs(fp.fft(hrw * input_signal)))[:len(k)]))
    high_res_win_amplitude_spectrum_ax.set_title("Амплитудный спектр сигнала с окном высокого разрешения")
    high_res_win_amplitude_spectrum_ax.set_xlabel("Частота, Гц")
    high_res_win_amplitude_spectrum_ax.set_ylabel("Амплитуда, Дб")

    # Вычисляем фазовый спектр сигнала с применение оконной функции низкого разрешения и строим график
    high_res_win_phase_spectrum = calculate_phase_spectrum(fp.fft(hrw * input_signal))
    high_res_win_phase_spectrum_ax = high_res_win_spectrum_window.add_subplot(212)
    high_res_win_phase_spectrum_ax.plot(frequency, high_res_win_phase_spectrum[:len(k)])
    high_res_win_phase_spectrum_ax.set_title("Фазовый спектр сигнала с окном высокого разрешения")
    high_res_win_phase_spectrum_ax.set_xlabel("Частота, Гц")
    high_res_win_phase_spectrum_ax.set_ylabel("Фаза, рад")

    # Добавляю пространство между графиками
    signal_window.tight_layout()
    fft_spectrum_window.tight_layout()
    low_res_win_spectrum_window.tight_layout()
    high_res_win_spectrum_window.tight_layout()

    signal_window.savefig(f"График {signal_name} дискретизированного сигнала.png")
    fft_spectrum_window.savefig(f"Спектры {signal_name} сигнала.png")
    low_res_win_spectrum_window.savefig(f"Спектры {signal_name} с применением оконной функции низкого разрешения.png")
    high_res_win_spectrum_window.savefig(f"Спектры {signal_name} с применением оконной функции высокого разрешения.png")

def signal(time):
    """ 
    Генерирует сигнал с заданными параметрами

    :time: отсчёты по х-оси
    :returns: сигнал с заданными параметрами
    """
    from numpy import sin, pi

    input_signal = A0 * sin(2 * pi * f0 * time) + A1 * sin(2 * pi * f1 * time) + A2 * sin(2 * pi * f2 * time)

    return input_signal 

def main():
    import numpy as np
    import matplotlib.pyplot as plt

    # Генерируем N отсчётов
    n = np.linspace(0, N - 1, N)

    # Генерируем дискретный сигнал
    sd = signal(n * delta_T) 

    # Расчитываем амплитудные и фазовые характеристики сигнала с окнами и без, и строим их графики
    signal_plot(sd, "обычного")

    # Генерируем шум с мат. ожиданием 0.03, среднеквадратическим отклонением 0.02 и длиной равной
    # длине дискретного сигнала
    noise = np.random.normal(0.03, 0.02, len(sd))

    # Повторяем расчёты и построения графика для сигнала с шумом
    signal_plot(sd + noise, "зашумлённого")

    # График шума
    plt.figure()
    plt.plot(range(0, len(sd), 1), noise)
    plt.xlabel("Отсчёты")
    plt.ylabel("Амплитуда, В")
    plt.title("Шум")
    plt.savefig("График шума")

    plt.show()

if __name__ == "__main__":
    main()
