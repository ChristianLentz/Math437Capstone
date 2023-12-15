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
        
        # set up evenly spaced, periodic data
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        xN = 9
        deltaX = (self.xmax - self.xmin) / xN 
        xvals = np.linspace(self.xmin, self.xmax, xN, endpoint = False)
        yvals = self.get_IC(xvals)

        # get a plot ready
        self.fig, self.ax = fig, ax

        # interpolate initial solution
        X_vals, Y_vals, max_val = self.transform(yvals)

        # plot the result
        self.ax.plot(X_vals, Y_vals.real, '-', color='blue')
        self.ax.set_ylim(0, max_val+1)
        
        # apply the wave equation 
        
        plt.show()
    
    def run(self):
        
        """
        Run the PDE solver and allow to plot/save! 
        """
    
    def get_IC(self, vals):
        
        """
        Generate an initial condition. 
        """
        
        rand_num = np.random.rand(1)
        if rand_num > 0.5:
            return np.exp(np.sin(vals))
        else:
            return np.exp(np.cos(vals))
    
    def transform(self, yvals): 
        
        """ 
        Perform a Fourier transform with scipy's fft(), then use a hardcoded
        version of the inverse fft to interpolate.  
        """
        
        # apply fft - with scipy
        cvals = fft(yvals) / yvals.size
        # apply inverse fft over a finer discretization of of xvals
        X_vals = np.linspace(self.xmin, self.xmax, 1000) 
        Y_vals = 0*X_vals+0j
        jvals = fftfreq(cvals.size, 1/cvals.size)
        max_val = 0
        for k in range(X_vals.size):
            new_val = (np.exp(1j*jvals*X_vals[k])*cvals).sum()
            if new_val.real > max_val: 
                max_val = new_val.real
            Y_vals[k] = new_val
        
        return X_vals, Y_vals, max_val
        
# ============================================================

def main():
    
    fig, ax = plt.subplots()
    b = 10*np.pi
    F = fourier([-b,b], [0,10], fig, ax)

if __name__ == "__main__": 
    main()      