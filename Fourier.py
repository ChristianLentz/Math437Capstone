""" 
This file contains code to approximate solutions to the wave equation 
in one dimension using the Fouier Spectral method via the fast Fourier 
Transform. 

Author: Christian Lentz
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq

class fourier: 
    
    def __init__(self, xbounds, tbounds, fig, ax): 
        
        """ 
        A PDE solver for the wave equation using Fourier spectral methods. 
        Checkout foureirExample.py to see a simpler version of this setup! 
        """
        
        # set up evenly spaced data - use the get IC to do this better! 
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        xN = 129
        deltaX = (self.xmax - self.xmin) / xN 
        xvals = np.linspace(self.xmin, self.xmax, xN, endpoint = False)
        yvals = np.random.rand(xN)

        # get a plot ready
        self.fig, self.ax = fig, ax
        self.ax.plot(xvals, yvals, 'o', color='red')

        # apply fft - with scipy
        cvals = fft(yvals) / yvals.size

        # apply inverse fft over a finer range of xvals
        X_vals = np.linspace(self.xmin, self.xmax, 1000) 
        Y_vals = 0*X_vals+0j
        jvals = fftfreq(cvals.size, 1/cvals.size)
        for k in range(X_vals.size):
            Y_vals[k] = (np.exp(1j*jvals*X_vals[k])*cvals).sum()

        # plot the result
        self.ax.plot(X_vals, Y_vals.real, '-', color='blue')
        self.ax.plot(X_vals, Y_vals.imag, '-', color='green')
        self.ax.set_ylim(-1,2)
        self.ax.set_title("Trig Interpolation with Fourier Transform")
        plt.show()
    
    def run(self):
        
        """
        Run the PDE solver and allow to plot/save! 
        """
    
    def get_IC(self):
        
        """
        Generate an initial condition. 
        """
    
    def transform(self, vals): 
        
        """ 
        Perform a Fourier transform with scipy's fft(), then use a hardcoded
        version of the inverse fft to interpolate.  
        """
        
# ============================================================

def main():
    
    fig, ax = plt.subplots()
    b = 10*np.pi
    F = fourier([-b,b], [0,10], fig, ax)

if __name__ == "__main__": 
    main()      