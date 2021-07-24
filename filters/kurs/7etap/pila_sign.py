def main():
    from scipy.signal import sawtooth
    import matplotlib.pyplot as plt
    import numpy as np

    waveform = 512 * (1 + sawtooth(2 * np.pi * np.linspace(0, 1, 40)))

    waveform = waveform.astype(np.int64)

    print((waveform))

    # plt.plot(waveform)
    # plt.show()


if __name__ == "__main__":
    main()
