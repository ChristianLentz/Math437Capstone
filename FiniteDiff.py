""" 
This file contains code to approximate solutions to the wave equation 
in two dimensions using finite differences. 

Author: Christian Lentz
"""

import numpy as np
import random as rand
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_toolkits.mplot3d import Axes3D 

class FiniteDiffs: 
    
    """ 
    A finite differences PDE solver which we use in FrontEnd.py
    """
    
    def __init__(self, xbounds, ybounds, tbounds, fig, ax):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # constants 
        self.c = 0.06
        
        # spatial and temporal bounds
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        self.ymin = ybounds[0]
        self.ymax = ybounds[1]
        self.tmin = tbounds[0]
        self.tmax = tbounds[1]
        
        # discretize spatial and temporal variables  
        xN = 75
        self.deltaX = (self.xmax - self.xmin) / xN
        self.xvals = np.linspace(self.xmin, self.xmax, xN) 
        yN = 75
        self.deltaY = (self.ymax - self.ymin) / yN
        self.yvals = np.linspace(self.ymin, self.ymax, yN) 
        tN = 50 
        self.deltaT = (self.tmax - self.tmin) / tN
    
        # initial conditions
        self.IC = np.array(self.get_IC())
        
        # array to hold previous and current solutions
        self.uprev = None             # u(x,y,T - dt)
        self.ucurrent = self.IC       # u(x,y,T)
        
        # counter for time steps
        self.tc = 0
        
        # figure to plot the surface of the wave
        self.fig = fig
        self.ax = ax
        
    def runTest(self): 
        
        """
        Called in the main method of this file. Will animate solutions to 
        a randomly generated in 
        """
        
        # prepare the grid and plot initial solution
        X, Y = np.meshgrid(self.xvals, self.yvals) 
        self.ax.plot_surface(X, Y, self.IC)
        self.ax.set_xlabel("x axis")
        self.ax.set_ylabel("y axis")
        self.ax.set_zlabel("z axis")
        
        # # plot successive solutions for each new time step
        # for tj in range(tN): 
        #     self.oneStep()
        #     self.ax.plot_surface(X, Y, self.ucurrent)
            
        # FINAL STEP IS TO ANIMATE!
        
        # animate the solution curves 
        # animation = ani.FuncAnimation(fig = self.fig, 
        #                                 func = FD.oneStep(),
        #                                 fargs = (img, ),  
        #                                 frames = 10, 
        #                                 interval = 50)
        
    def get_dt(self):  
               
        """
        Use the CLT condition to find a value for dt given dx and dy. 
        """
        
        # this may not be included, we'll see
    
    def get_IC(self): 
        
        """
        Randomly generate an initial condition. 
        We use random linear combintions of normal distributions. 
        """
        
        # generate random x and y paths as initial conditions
        ICx = self.getRandHillyPath(self.xvals, self.xmax, self.xmin)
        ICy = self.getRandHillyPath(self.yvals, self.ymax, self.ymin)
        
        # combine ICx and ICy into a surface
        IC = []
        for xval in ICx:
            curr_x_arr = []
            for yval in ICy: 
                curr_x_arr.append(xval+yval)
            IC.append(curr_x_arr)
            
        return IC
    
    def getRandHillyPath(self, vals, maximum, minimum):
        
        """ 
        Generate a random hilly path. We use this to generate an initial 
        condition for y and x, and then combine these to get a random 
        surface as our initial wave. 
        """
        
        ic = np.zeros(vals.shape)
        for hill in range(rand.randrange(5)):
            # random spread for each hill  
            sigma = np.random.random() * (maximum/2 - minimum)
            # random center for each hill 
            mu = np.random.random() * (maximum - minimum)
            # add new hill to the initial x condition
            ic = ic + self.getRandGaussian(vals, sigma, mu)
            
        return ic
    
    def getRandGaussian(self, vals, sigma, mu):
        
        """
        Generate a Guassian (normal) distribution with random spread and center. 
        """
        return np.array((1 / (np.sqrt(2 * np.pi * sigma**2))) * np.exp((-1/2) * ((vals-mu) / sigma)**2))
        
    def oneStep(self): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - xl = x - dx = the previous x coordinate
            - xc = x = the current x coordinate
            - xr = x + dx = the next x coordinate 
            
            - yl = y - dy = the previous y coordinate 
            - yc = y = the current y coordinate 
            - yr = y + dy = the next y coordinate
            
            - tp = t - dt = the previous time step
        """

        unew = 0 * self.ucurrent
        for j in range(unew.size):
            
            # get tp if current time step is not the initial 
            if self.tc != 0: 
                tp = self.uprev[j]
            else: 
                tp = 0
                
            # left boundary -- UPDATE THIS TO BE PERIODIC
            if j == 0:
                unew[j] = 0
                
            # right boundary -- UPDATE THIS TO BE PERIODIC
            elif j == unew.size-1:
                # get components of difference quotient 
                xl = self.ucurrent[j-1]
                xc = self.ucurrent[j]
                # compute difference quotient
                unew[j] = 2*xc - tp
                unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(0 + xl - 2*xc)
            
            else: 
                # get components of difference quotient 
                xl = self.ucurrent[j-1]
                xc = self.ucurrent[j]
                xr = self.ucurrent[j+1]
                # compute difference quotient
                unew[j] = 2*xc - tp 
                unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(xr + xl - 2*xc) 
        
        # update solutions and time step
        self.uprev = self.ucurrent
        self.ucurrent = unew
        self.tc+=1

# ==============================================================

def main():
    
    """
    The main method here is for testing the finite differences code. 
    """ 
    
    # instantiate finite diff solver 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    FD = FiniteDiffs([0,75], [0, 75], [0,25], fig, ax)
    
    FD.runTest()
        
    plt.show()

if __name__ == '__main__':
    main()