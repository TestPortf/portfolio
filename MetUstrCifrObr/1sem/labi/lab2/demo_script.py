import scipy.signal as sig
import matplotlib.pyplot as plt

win = sig.windows.nuttall(1024)
plt.plot(range(0, 1024, 1), win)
plt.show()

