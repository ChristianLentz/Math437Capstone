"""
A "hello world" of the Fourier Spectral method. 

This just genertes points and does some trig interpolation using using 
the Fourier transform!  
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq

# set up evenly spaced data 
xmin = 0
xmax = 2*np.pi
xN = 9
deltaX = (xmax - xmin) / xN 
xvals = np.linspace(xmin, xmax, xN, endpoint = False)
yvals = np.random.rand(xN)

# get a plot ready
fig, ax = plt.subplots() 
ax.plot(xvals, yvals, 'o', color='red')

# apply fft 
cvals = fft(yvals) / yvals.size

# apply inverse fft over a finer range of xvals
X_vals = np.linspace(0, 2*np.pi, 1000) 
Y_vals = 0*X_vals+0j
jvals = fftfreq(cvals.size, 1/cvals.size)
for k in range(X_vals.size):
    Y_vals[k] = (np.exp(1j*jvals*X_vals[k])*cvals).sum()

# plot the result
ax.plot(X_vals, Y_vals.real, '-', color='blue')
ax.plot(X_vals, Y_vals.imag, '-', color='green')
ax.set_ylim(-1,2)
ax.set_title("Trig Interpolation with Fourier Transform")
plt.show()