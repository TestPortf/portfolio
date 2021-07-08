import numpy as np
import sounddevice as sd
import time

# Samples per second
sps = 44100

# Frequency / pitch
freq_hz = 440.0

# Duration
duration_s = 5.0

# Attenuation so the sound is reasonable
atten = 1 

# NumpPy magic to calculate the waveform
waveform = [0, 231, 464, 697, 929, 1162, 1395, 1627, 1859, 2001, 1768, 1536, 1303, 1070, 839, 606, 373, 141, -91, -324, -556, -788, -1021, -1254, -1486, -1719, -1952, -1911, -1679, -1446, -1213, -981, -748, -515, -283]
waveform_quiet = waveform * atten

# Play the waveform out the speakers
sd.play(waveform_quiet, sps)
time.sleep(duration_s)
sd.stop()