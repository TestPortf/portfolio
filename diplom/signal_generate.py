dac_bit = 2**12 - 1
discretization_freq = 22000
discretization_period = 1/discretization_freq
phase_accamulator_bitness = 32
phase_accamulator = 2**phase_accamulator_bitness
wavetable_length_bitness = 8
wavetable_length = 2**wavetable_length_bitness
accumulator_comparison_shift = phase_accamulator_bitness - wavetable_length_bitness


def print_c_array(generated_wavetable, name, columns=16):
    print(name + ' {' + str(generated_wavetable[0]), end=', ')
    for current_signal_count in range(1, len(generated_wavetable) - 1):
        if (current_signal_count % columns) == 0:
            print('\\')
            print(' ' * len(name), end='  ')
        if generated_wavetable[current_signal_count] < 0:
            if abs(generated_wavetable[current_signal_count]) < 10:
                print('  ', generated_wavetable[current_signal_count], end=', ')
            elif abs(generated_wavetable[current_signal_count]) < 100:
                print(' ', generated_wavetable[current_signal_count], end=', ')
            elif abs(generated_wavetable[current_signal_count]) < 1000:
                print('', generated_wavetable[current_signal_count], end=', ')
            else: print(generated_wavetable[current_signal_count], end=', ')
        else:
            if generated_wavetable[current_signal_count] < 10:
                print('   ', generated_wavetable[current_signal_count], end=', ')
            elif generated_wavetable[current_signal_count] < 100:
                print('  ', generated_wavetable[current_signal_count], end=', ')
            elif generated_wavetable[current_signal_count] < 1000:
                print(' ', generated_wavetable[current_signal_count], end=', ')
            else: print('', generated_wavetable[current_signal_count], end=', ')
    print(generated_wavetable[-1], end='};')
    print('')


def meandr_gen(t, freq):
    import numpy as np
    meandr_arr = np.sign(np.sin(2 * np.pi * freq * t))
    #meandr_arr[meandr_arr < 0] = 0
    return dac_bit/2 * meandr_arr

def triangle_wave(t, freq):
    import numpy as np
    return ((2 * dac_bit)/np.pi * np.arcsin(np.sin(((2 * np.pi) * t *freq))) )/2


def sin_wave(t, freq):
    import numpy as np
    return (dac_bit * np.sin(2 * np.pi * t * freq) )/2 


def sawtooth_wave(t, freq):
    import numpy as np
    return ((2 * dac_bit)/np.pi * np.arctan(np.tan(((2 * np.pi) * t *freq/2))))/2


def amplitude_spectrum_calculate(signal_to_analisis):
    import numpy as np
    import matplotlib.pyplot as plt

    number_of_signal_repetitions = 80
    repeated_signal = np.tile(signal_to_analisis, number_of_signal_repetitions)

    fft_signal = np.fft.fft(repeated_signal)
    frequency = np.linspace(0, len(fft_signal), len(fft_signal))
    frequency = frequency/(len(fft_signal) * discretization_period)
    plt.plot(frequency[:int(len(repeated_signal)/2)], abs(fft_signal)[:int(len(repeated_signal)/2)])
    plt.xlabel('Частота, Гц')
    plt.ylabel('Амплитуда')


def get_phase_increment(expected_frequency):
    return int(phase_accamulator * expected_frequency/discretization_freq)


def signal_interpolation(current_phase_index, wavetable):
    comparison_mask = (1 << accumulator_comparison_shift) - 1

    nearest_table_index = current_phase_index >> accumulator_comparison_shift
    y_0 = wavetable[nearest_table_index]
    y_1 = wavetable[(nearest_table_index + 1) % wavetable_length]

    vect = phase_accamulator/wavetable_length
    return (y_0 * (vect - comparison_mask) + y_1 * comparison_mask)/vect


def signal_interpolation_without_div(current_phase_index, wavetable):
    comparison_mask = (1 << accumulator_comparison_shift) - 1

    nearest_table_index = current_phase_index >> accumulator_comparison_shift
    y_0 = wavetable[nearest_table_index]
    y_1 = wavetable[(nearest_table_index + 1) % wavetable_length]

    phase_accum_and_table_dif = (comparison_mask&current_phase_index) 
    return (y_0 * ((1 << accumulator_comparison_shift) - phase_accum_and_table_dif) + y_1 * phase_accum_and_table_dif) >> accumulator_comparison_shift


def calculate_phase_increment_for_octave(octave_frequences, octave_name):
    phase_array = []
    for note in range(0, len(octave_frequences)):
        phase_array.append(get_phase_increment(octave_frequences[note]))
    print_c_array(phase_array, f'unsigned long int phase_incriment_{octave_name}_oct[13] = ')


def main():
    import numpy as np
    import matplotlib.pyplot as plt
    import sys

    small_octave_frequences = [130.813, 138.591, 146.832, 155.563, 164.814, 174.614, 184.997, 195.998, 207.652, 220, 233.082, 246.942, 261.626]
    first_octave_frequences = [261.626, 277.183, 293.655, 311.127, 329.628, 349.228, 369.994, 391.995, 415.305, 440, 466.164, 493.883, 523.251]
    second_octave_frequences= [523.251, 554.365, 587.330, 622.254, 659.255, 698.456, 739.989, 783.991, 830.609, 880, 932.328, 987.767, 1046.50]
    third_octave_frequences = [1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760,1864.66, 1975.53, 2093]

    time = np.linspace(0, wavetable_length * discretization_period, wavetable_length)
    sample = sin_wave(time, discretization_freq/wavetable_length).astype(int)
    print_c_array(sample, 'int sin_wave[256] = ')
    sample = meandr_gen(time, discretization_freq/wavetable_length).astype(int)
    print_c_array(sample, 'int meandr_wave[256] = ')
    sample = sawtooth_wave(time, discretization_freq/wavetable_length).astype(int)
    print_c_array(sample, 'int sawtooth_wave[256] = ')
    sample = triangle_wave(time, discretization_freq/wavetable_length).astype(int)
    print_c_array(sample, 'int triangle_wave[256] = ')

    note_frequency = second_octave_frequences[3]
    ph_m = get_phase_increment(note_frequency)
    ph_i = 0
    probe_number = int(phase_accamulator/ph_m)
    different_frequency_signal = []
    for current_signal_count in range(0, probe_number):
        different_frequency_signal.append(signal_interpolation_without_div(ph_i, sample))
        ph_i = (ph_i + ph_m) % phase_accamulator 

    print(different_frequency_signal)

    # plt.figure()
    # plt.plot(sample)
    # plt.plot(different_frequency_signal[:int(phase_accamulator/ph_m)])
    # plt.legend(['Табличный сигнал', 'Сигнал полученный из табличного методом интерполяции'])
    # plt.title(f'Выбраная частота передисретизации сигнала = {note_frequency}')

    # plt.figure()
    # amplitude_spectrum_calculate(16*sample)
    # amplitude_spectrum_calculate(different_frequency_signal)
    # plt.legend(['Табличный сигнал', 'Сигнал полученный из табличного методом интерполяции'])
    # plt.title(f'Выбраная частота передисретизации сигнала = {note_frequency}')

    # calculate_phase_increment_for_octave(small_octave_frequences, 'small')
    # calculate_phase_increment_for_octave(first_octave_frequences, 'first')
    # calculate_phase_increment_for_octave(second_octave_frequences, 'second')
    # calculate_phase_increment_for_octave(third_octave_frequences, 'third')


    plt.show()

if __name__ == "__main__":
    main()