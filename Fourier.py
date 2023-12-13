# """ 
# This file contains code to approximate solutions to the wave equation 
# in one dimension using the Fouier Spectral method via the fast Fourier 
# Transform. 

# Author: Christian Lentz
# """

# class fourier: 
    
#     def __init__(self): 
        
#         """ 
#         A PDE solver for the wave equation using Fourier spectral methods. 
#         """
        
#         self.xmin = 0
#         self.xmax = 2*pi
#         xN = 2 * (self.xmax - self.xmin)
#         self.deltaX = (self.xmax - self.xmin) / 8 
#         self.xvals = np.linspace(self.xmin, self.xmax, xN, endpoint = False) 
#         tN = 1000 * (self.tmax - self.tmin) 
#         self.deltaT = (self.tmax - self.tmin) / tN
#         self.XC = ((self.c**2)*(self.deltaT**2))/(self.deltaX**2)
        
#         print("I worked")
        
    
#     def transform(self, vals): 
        
#         """ 
#         Perform a Fourier transform on some vecotr using the Fourier matrix. 
#         """
        
#     def main():
        
#         fourier = fourier()

# if __name__ == "__main__": 
#     main()

import numpy as np
import math as math 
import random as rand
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq

# set up evenly spaced data
xmin = 0
xmax = 2*np.pi
xN = 9
deltaX = (xmax - xmin) / 9 
xvals = np.linspace(xmin, xmax, xN, endpoint = False)
yvals = np.random.rand(xN)
# get a plot ready
fig, ax = plt.subplots() 
ax.plot(xvals, yvals, 'o', color='red')
# apply fft 
cvals = fft(yvals) / yvals.size
print(cvals)
# apply inverse fft to get data to plot 
# real part should be small!
X_vals = np.linspace(0, 2*np.pi, 1000) 
Y_vals = 0*X_vals+0j
jvals = fftfreq(cvals.size, 1/cvals.size)
print(jvals)
# jvals = np.arange(8)
for k in range(X_vals.size):
    Y_vals[k] = (np.exp(1j*jvals*X_vals[k])*cvals).sum()
ax.plot(X_vals, Y_vals.real, '-', color='blue')
ax.plot(X_vals, Y_vals.imag, '-', color='green')
ax.set_ylim(-1,2)
plt.show()

# tN = 1000 * (tmax - tmin) 
# deltaT = (tmax - tmin) / tN
# XC = ((c**2)*(deltaT**2))/(deltaX**2)


        