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
    
    def __init__(self, xbounds, tbounds, fig, ax, ybounds = None):
        
        """ 
        Constructor for the finite differences solver. 
        """
        
        # constants and variables  
        self.c = 1.5
        if ybounds == None:
            self.is1D = True 
        else: 
            self.is1D = False
        
        # spatial and temporal bounds
        self.xmin = xbounds[0]
        self.xmax = xbounds[1]
        self.tmin = tbounds[0]
        self.tmax = tbounds[1]
        
        # discretize spatial and temporal variables  
        xN = 75
        self.deltaX = (self.xmax - self.xmin) / xN
        self.xvals = np.linspace(self.xmin, self.xmax, xN) 
        tN = 2 * (self.tmax - self.tmin) 
        self.deltaT = (self.tmax - self.tmin) / tN
        
        # add second spatial dimension if necessary
        if not self.is1D:
            self.ymin = ybounds[0]
            self.ymax = ybounds[1]
            yN = 75
            self.deltaY = (self.ymax - self.ymin) / yN
            self.yvals = np.linspace(self.ymin, self.ymax, yN) 
        
        # initial conditions and previous solutions
        # xprev and yprev together give u(x,y,t - dt)
        if self.is1D: 
            self.ICx = self.get_IC(self.is1D)
        else:
            self.ICx, self.ICy = self.get_IC(self.is1D)
            self.yprev = None
            self.ycurr = self.ICy
            self.IC = self.projectToSurface(self.ICx, self.ICy)
        self.xprev = None
        self.xcurr = self.ICx
       
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
        
        # plot 1D solutions 
        if self.is1D:
            self.plotSoln1D(self.ICx)
            for tj in range(self.tmax - 1):
                self.oneStep1D()
                self.plotSoln1D(self.xcurr)
        # plot 2D solutions 
        else:
            self.plotSoln2D(self.IC)
            for tj in range(self.tmax-1):
                self.oneStep2D()
                self.plotSoln2D(self.projectToSurface(self.xcurr, self.ycurr))
                
    def get_dt(self):  
               
        """
        Use the CLT condition to find a value for dt given dx and dy. 
        """
        
        # this may not be included, we'll see
    
    def get_IC(self, is1D): 
        
        """
        Randomly generate an initial condition. 
        We use random linear combintions of normal distributions. 
        """
        
        # generate random x and/or y paths as initial conditions
        if is1D: 
            return self.getRandHillyPath(self.xvals, self.xmax, self.xmin)
        else: 
            ICx = self.getRandHillyPath(self.xvals, self.xmax, self.xmin)
            ICy = self.getRandHillyPath(self.yvals, self.ymax, self.ymin)
            return ICx, ICy
    
    def getRandHillyPath(self, vals, maximum, minimum):
        
        """ 
        Generate a random hilly path. We use this to generate an initial 
        condition for y and x, and then combine these to get a random 
        surface as our initial wave. 
        """
        
        ic = np.zeros(vals.shape)
        for hill in range(rand.randrange(10)):
            # random spread for each hill  
            sigma = (1 - np.random.random()) * (maximum/2 - minimum)
            # random center for each hill 
            mu = (1 - np.random.random()) * (maximum - minimum)
            # add new hill to the initial x condition
            ic = ic + self.getRandGaussian(vals, sigma, mu)
            
        return ic
    
    def getRandGaussian(self, vals, sigma, mu):
        
        """
        Generate a Guassian (normal) distribution with random spread and center. 
        """
        return np.array((1 / (np.sqrt(2 * np.pi * sigma**2))) * np.exp((-1/2) * ((vals-mu) / sigma)**2))
    
    def projectToSurface(self, xvals, yvals): 
        
        """
        From two sequences which represent points along the x and y axis, project 
        these into a surface. 
        """
        
        surf = []
        for xval in xvals:
            curr_x_arr = []
            for yval in yvals: 
                curr_x_arr.append(xval+yval)
            surf.append(curr_x_arr)
            
        return np.array(surf)
        
    def oneStep1D(self): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - xl = x - dx = the previous x coordinate
            - xc = x = the current x coordinate
            - xr = x + dx = the next x coordinate 
            
            - tp = t - dt = the previous time step
        """

        unew = 0 * self.xcurr
        
        for j in range(unew.size):
            
            # get tp if current time step is not the initial 
            if self.tc != 0: 
                tp = self.xprev[j]
            else: 
                tp = 0
                
            # left boundary -- periodic boundary condition
            if j == 0:
                unew[j] = self.xcurr[self.xcurr.size-1]
                
            # right boundary -- periodic boundary condition
            elif j == unew.size-1:
                unew[j] = self.xcurr[0]
            else: 
                # get components of difference quotient 
                xl = self.xcurr[j-1]
                xc = self.xcurr[j]
                xr = self.xcurr[j+1]
                # compute difference quotient
                unew[j] = 2*xc - tp 
                unew[j] += (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(xr + xl - 2*xc) 
        
        # update solutions and time step
        self.xprev = self.xcurr
        self.xcurr = unew
        self.tc+=1
        
    def plotSoln1D(self, vals): 
        
        """
        Plot a 1D solution to the wave equation. 
        """
        
        self.ax.plot(self.xvals, vals, '-')
        self.ax.set_xlabel("x axis")
        self.ax.set_ylabel("y axis")
        
    def oneStep2D(self): 
        
        """ 
        Run one step of the PDE solver by computing a solution for a single time.
        Each piece of the difference quotient that we need is written as follows: 
        
            - xl = x - dx = the previous x coordinate
            - xc = x = the current x coordinate
            - xr = x + dx = the next x coordinate 
            
            - yl = y - dy = the previous y coordinate 
            - yc = y = the current y coordinate 
            - yr = y + dy = the next y coordinate
            
            - tpx = t - dt = the previous time step for x
            - tpy = t - dt = the previous time step for y
        """

        xnew = 0 * self.xcurr
        ynew = 0 * self.ycurr
         
        for j in range(self.xcurr.size):
            
            # get tp if current time step is not the initial 
            if self.tc != 0: 
                tpx = self.xprev[j]
                tpy = self.xprev[j]
            else: 
                tpx = 0
                tpy = 0
                
            # left boundary -- periodic boundary condition
            if j == 0:
                xnew[j] = self.xcurr[self.xcurr.size-1]
                ynew[j] = self.ycurr[self.ycurr.size-1]
                
            # right boundary -- periodic boundary condition 
            elif j == xnew.size-1:
                xnew[j] = self.xcurr[0]
                ynew[j] = self.ycurr[0]
            
            else: 
                # get components of difference quotient 
                xl = self.xcurr[j-1]
                xc = self.xcurr[j]
                xr = self.xcurr[j+1]
                yl = self.ycurr[j-1]
                yc = self.ycurr[j]
                yr = self.xcurr[j+1]
                # compute difference quotient
                xnew[j] = 2*xc - tpx + (((self.c**2)*(self.deltaT**2))/(self.deltaX**2))*(xr + xl - 2*xc) 
                ynew[j] = 2*yc - tpy + (((self.c**2)*(self.deltaT**2))/(self.deltaY**2))*(yr + yl - 2*yc)
        
        # update solutions and time step
        self.xprev = self.xcurr
        self.yprev = self.ycurr
        self.xcurr = xnew
        self.ycurr = ynew
        self.tc+=1
        
    def plotSoln2D(self, vals): 
        
        """
        Plot a 2D solution to the wave equation. 
        """
        
        X, Y = np.meshgrid(self.xvals, self.yvals) 
        self.ax.plot_surface(X, Y, vals)
        self.ax.set_xlabel("x axis")
        self.ax.set_ylabel("y axis")
        self.ax.set_zlabel("z axis")

# ==============================================================

def main():
    
    """
    The main method here is for testing the finite differences code. 
    """ 
    
    # instantiate finite diff solver 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    FD = FiniteDiffs([0,75], [0,100], fig, ax, ybounds = [0, 75])
    
    # ybounds = [0, 75]
    FD.runTest()
    plt.show()

if __name__ == '__main__':
    main()